module SingleHole

using OrderedCollections
using SparseArrays
using KrylovKit
using JSON

"""
`struct System` immutable structure for system input parameters:
# Fields
*   `size::Int64`: number of lattice sites
*   `momentum::Int64`: index indicating which momentum sector should be calculated. Index `momentum` corresponds to momentum `k = 2π * momentum / size`.
*   `magnetization::Int64`: index in range `0:floor(size/2)` indicating which magnetization sector should be calculated. Magnetization is taken without sign.
*   `coupling::Float64`: value of the coupling constant J. For ferromagnet `coupling < 0`, while for antiferromagnet `coupling > 0`.
*   `interaction::FLoat64`: parameter for scaling magnon-magnon interactions. For pure Heisenberg model `interaction = 1.0`.
"""
struct System
    size::Int64
    momentum::Int64
    magnetization::Int64
    coupling::Float64
    interaction::Float64
    intrachain::Float64
    isotropy::Float64
    hopping::Float64
end

"`Basis === OrderedDict{Int64, Int64}`"
Basis = OrderedDict{Int64, Int64}

"""
`mutable struct LinearCombination`: structure for storing result of operators action on states belonging to `basis::Basis`
# Fields
*   `state::Int64`: magnon configuration in binary representation written as decimal number
*   `coefficient::Vector{Complex{Float64}}`: coeffcient multiplying state in the linear combination
"""
mutable struct LinearCombination
    state::Vector{Int64}
    coefficient::Vector{Complex{Float64}}
end

"`Model === SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

"""
    run()

Run t-J model diagonalization procedure.
"""
function run(input::Union{Missing, OrderedDict} = missing; howmany = 10, factor = true)
    system::System = if input === missing
        readInput()
    else
        System(
            input["system size"],
            input["momentum sector"],
            input["magnetization sector"],
            input["coupling constant"],
            input["magnon interaction"],
            input["intrachain coupling"],
            input["intrachain isotropy"],
            input["hopping constant"]
        )
    end
    println()
    @time basis::Basis = makeBasis(system)
    @time model::Model = makeModel(basis, system)
    if factor
        @time factorization = factorize(model, howmany = howmany, kryldim = max(30, 3*howmany)) # put 10 for DOS
        return system, basis, model, factorization
    else
        return system, basis, model
    end
end

"Read `input.json` file and return `System` structure with input data. It requires `input.json` file to be located in the current working directory."
function readInput()::System
    path = "./single_hole/code/julia/"
    input = JSON.parsefile(path * "input.json", use_mmap = false) # use_mmap = false is a workaroud to ensure one can change JSON file without restarting Julia
    return System(
        input["system size"],
        input["momentum sector"],
        input["magnetization sector"],
        input["coupling constant"],
        input["magnon interaction"],
        input["intrachain coupling"],
        input["intrachain isotropy"],
        input["hopping constant"]
    )
end

"""
    makeBasis(system::System) -> Basis

Return `Basis === Dict{Int64, Int64}` dictionary representing basis of given magnetization and momentum sector specified by `system.magnetization` and `system.momentum` respectively. Each index in `Basis` corresponds to state value, and each value in `Basis` corresponds to position in the basis.
"""
function makeBasis(system::System)::Basis
    if isodd(system.size)
        error("Requested 'system.size' is odd! Only even sizes are supported.")
    end
    if system.magnetization < 0 || 2 * system.magnetization > system.size ^ 2
        error("Wrong magnetization sector in the input file! Possible subspaces are denoted by integers from 0 to N/2.")
    end

    ### Note: in the calculated basis representative states
    ###       have hole placed at site = system.size^2

    ### number of spins pointing up
    nSpinsUp::Int64 = div(system.size ^ 2, 2) - system.magnetization
    ### note: in the code `spin up === 1`, `spin down === 0`

    ### calculate magnetic subspace size
    subspaceSize::Int64 = binomial(system.size ^ 2, nSpinsUp)

    ### get first state (i.e. with lowest index in binary basis)
    ### note: `1 << n == 2^n`, but former is faster
    state::Int64 = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))

    ### set sublattice rotation masks
    mask = getMask(system)

    ### initialize the basis
    basis::Basis = Basis()

    index = 0
    ### iterate over states within given magnetization subspace
    for _ in 1:subspaceSize
        ### perform sublattice rotation
        magnonState = sublatticeRotation(state, mask)
        ### now 1 represents magnon, 0 stands for empty site

        ### add magnon state to basis
        push!(basis, magnonState => (index += 1))
        ### get next state
        state = getNextState(state, system)
    end

    ### repeat for reversed magnetization (if exists)
    if system.magnetization != 0
        ### number of spins pointing up
        nSpinsUp = div(system.size ^ 2, 2) + system.magnetization
        ### note: in the code `spin up === 1`, `spin down === 0`

        ### get first state (i.e. with lowest index in binary basis)
        ### note: `1 << n == 2^n`, but former is faster
        state = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))

        ### iterate over states within given magnetization subspace
        for _ in 1:subspaceSize
            ### perform sublattice rotation
            magnonState = sublatticeRotation(state, mask)
            ### now 1 represents magnon, 0 stands for empty site

            ### add magnon state to basis
            push!(basis, magnonState => (index += 1))
            ### get next state
            state = getNextState(state, system)
        end
    end

    return basis
end

"""
    getMask(system::System) -> Int64

Return 2D AF mask.
"""
function getMask(system::System)::Int64
    mask = 0
    for i = 0:(system.size - 1)
        chainMask = sum(1 << k for k in 0:2:(system.size - 1))
        if mod(i, 2) != 0
            chainMask = 2^system.size - 1 - chainMask
        end
        mask += (1 << (system.size * i)) * chainMask
    end
    return mask
end

"""
    sublatticeRotation(state::Int64, mask::Int64) -> Int64

Reverse bits according to mask.
"""
function sublatticeRotation(state::Int64, mask::Int64)::Int64
    return xor(state, mask)
end

"""
    getNextState(system::System, state::Int64) -> Int64

Return `newState::Int64` next to `state::Int64` within given magnetization sector specified by `system.magnetization`.
Two states, `state` and `newState`, are next to each other withing given magentization sector if and only if there is no state that has value between `state` and `newState`.
"""
function getNextState(state::Int64, system::System)::Int64
    count = 0
    ### loop over <i,j> site pairs (bonds)
    for i in 0:(system.size^2 - 2)
        j = i + 1

        ### get numeric value at i and j bit positions
        iValue, jValue = 1 << i, 1 << j

        ### check if there is bit 1 at site i
        if state & iValue > 0
            ### check if there is bit 1 at site j
            if state & jValue > 0
                ### clear bit at site i
                state &= ~iValue

                ### raise counter of cleared bits
                count += 1
            else
                ### swap bits at i and j
                state = xor(state, iValue + jValue)

                ### add cleared bits to lowest positions
                state += count == 0 ? 0 : sum(n -> 1 << n, 0:(count-1))

                ### return new state
                return state
            end
        end
    end
    ### if loop ends without returning new state we return -1
    ### as there is no next state in the requestes subspace
    return -1
end

"""
    bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) -> Int

Cyclic bit shift for calculationg bit translations with periodic boundary conditions.

# Arguments
*   `s::Int` - value which binary representation will be shifted
*   `l::Int` - size of the cycle (total number of bits in the cycle)
*   `f::Bool` - `true ->` move forward (~mult by `2`); `false ->` backward (~div by `2`); ~mult and ~div mean multiplication and division with preservance of periodic boundary conditions within cycle size `l`
*   `hb::Int` [optional] - highest bit (for speed up put the value of 2^(l-1))
*   `hv::Int` [optional] - highest value (for speed put (2^l)-1)
"""
@inline bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) = f ? 2s - div(s, hb) * hv : div(s, 2) + rem(s, 2) * hb

@inline function bitmov2(s::Int, d::Int, f::Bool = false; hb::Int = 1 << (d - 1), hv::Int = (1 << d) - 1)::Int64
    result = 0
    for i in 0:(d-1)
        take = (2^d)^(i+1) - (2^d)^i
        val = (s & take) >> (d * i)
        result += bitmov(val, d, f, hb = hb, hv = hv) << (d * i)
    end
    return result
end

function showState(state::Int64, system::System)
    result = Array{Int64, 2}(undef, system.size, system.size)
    for i in 0:(system.size-1)
        d = digits(state, base = 2, pad = system.size * system.size)
        result[i+1, :] = d[(1 + i * system.size):(system.size + i * system.size)]
    end
    return result
end

"""
    act(operator::Function, state::Int64, basis::Basis, system::System) -> LinearCombination

Apply `operator` to `state` belonging to `basis` and returns `LinearCombination  === Dict{Int64, Complex{Float64}}` representing states with their coefficients.
"""
function act(operator::Function, state::Int64, basis::Basis, system::System)::LinearCombination
    return operator(state, basis, system)
end

"""
    hamiltonian(state::Int64, basis::Basis, system::System) -> LinearCombination

Apply Hamiltonian to `state` written in Sz momentum `basis` obtained for input `system` parameters. Returns `LinearCombination  === Dict{Int64, Complex{Float64}}` representing resulting states with their coefficients.
"""
function hamiltonian(state::Int64, basis::Basis, system::System)::LinearCombination
    ### initialize result as linear combination
    result = LinearCombination(fill(state, 2 * system.size^2 + 1), zeros(Complex{Float64}, 2 * system.size^2 + 1))

    ### check if initial state belongs to basis
    if haskey(basis, state)
        ### initialize some constants for faster evaluation
        l::Int = system.size
        highestBit::Int = 1 << (l - 1)
        highestValue::Int = (1 << l) - 1

        ### initialize ik for faster exponent calculations
        ik::Complex{Float64} = 2.0 * pi * im * system.momentum / system.size

        ### state periodicity
        periodicity = system.size

        ### loop over lattice sites without hole
        ### case 1: interchain
        for i in 1:(system.size^2 - 2)
            j = mod1(i + 1, system.size) + system.size * div(i - 1, system.size)

            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))

            ## work out off-diagonal coeffcients
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            if iBit == jBit
                ### if two neighbouring spins are different then flip those spins
                ### here we have magnon language so if 2 sites are both empty (or
                ### occupied) by magnons then create (or annihilate) magnons
                newState = xor(state, iValue + jValue)

                ### calculate matrix coefficient
                coefficient = 0.5   # there is no exponantial since distance to representative is 0

                ### put new state in linear combination
                ### and set its corresponding coefficient
                result.state[i + 1] = newState
                result.coefficient[i + 1] = coefficient * system.coupling
            end

            #### diagonal coefficient
            ## comment: after rotation bits represent magnons (0 -> no magnon, 1 -> magnon present)
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            result.coefficient[1] -= (0.25 - 0.5 * (iBit + jBit) + system.interaction * iBit * jBit) * system.coupling
        end

        ### n_i n_j terms
        result.coefficient[1] += 0.5 * system.coupling

        ### workout sites/bonds around the hole for off-diagonal coeffcients
        i = system.size^2 # hole position
        for j in [system.size^2 - 1, system.size^2 - system.size + 1]
            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)

            newState = state
            ## if no magnon at site of arrival
            if iBit == jBit
                newState = xor(state, iValue + jValue)
            end
            ### calculate matrix coefficient
            coefficient = 1.0   # init
            if i - j > 1    # i - j > 1 -> j === system.size^2 - system.size + 1 (we jump forward)
                ### so we translate backward
                ### and exp is taken with positive coeffcient exp(+)
                newState = bitmov2(newState, l, false, hb = highestBit, hv = highestValue)
                coefficient *= exp(ik)
            else
                ### else hole is at system.size-1 position (we jumped backward)
                ### so we translate forward and we have exp(-)
                newState = bitmov2(newState, l, true, hb = highestBit, hv = highestValue)
                coefficient *= exp(-ik)
            end

            ### create a new entry in linear combination
            ### and set its corresponding coeffcient
            j = mod(j, system.size - 1) # trick, but working :]
            result.state[system.size^2 + j] = newState
            result.coefficient[system.size^2 + j] = coefficient * system.hopping
        end

        ### case 1: intrachain
        for i in 1:(system.size^2 - 1)
            j = mod1(i + system.size, system.size^2)

            if j == system.size^2
                continue
            end

            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))

            ## work out off-diagonal coeffcients
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            if iBit == jBit
                ### if two neighbouring spins are different then flip those spins
                ### here we have magnon language so if 2 sites are both empty (or
                ### occupied) by magnons then create (or annihilate) magnons
                newState = xor(state, iValue + jValue)

                ### calculate matrix coefficient
                coefficient = 0.5  # there is no exponantial since distance to representative is 0

                ### put new state in linear combination
                ### and set its corresponding coefficient
                result.state[i + system.size^2 + 1] = newState
                result.coefficient[i + system.size^2 + 1] = coefficient * system.intrachain * system.isotropy
            end

            #### diagonal coefficient
            ## comment: after rotation bits represent magnons (0 -> no magnon, 1 -> magnon present)
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            result.coefficient[1] -= (0.25 - 0.5 * (iBit + jBit) + system.interaction * iBit * jBit) * system.intrachain
        end

        ### n_i n_j terms
        result.coefficient[1] += 0.5 * system.intrachain
    end

    ### return resulting linear combination
    return result
end

"""
    makeModel(basis::Basis, system::System) -> Model

Calculate dense matrix of the `model` Hamiltonian. Returns `Model === Array{Complex{Float64},2}`.
"""
function makeModel(basis::Basis, system::System)::Model
    subspaceSize = length(basis)
    linearCombinationLength = 2 * system.size^2 + 1
    I = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    J = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    V = Vector{Complex{Float64}}(undef, linearCombinationLength * subspaceSize)
    for (state, index) in basis
        linearCombination::LinearCombination = act(hamiltonian, state, basis, system)
        for it in 1:linearCombinationLength
            I[(index - 1) * linearCombinationLength + it] = index
            J[(index - 1) * linearCombinationLength + it] = basis[linearCombination.state[it]]
            V[(index - 1) * linearCombinationLength + it] = linearCombination.coefficient[it]
        end
    end
    return dropzeros!(sparse(I, J, V, subspaceSize, subspaceSize, +))
end

"""
    factorize(model::Model, [howmany = 1, which = :SR])

Compute eigenvalues (by default with smallest real part) and their corresponding eigenvectors.
"""
function factorize(model::Model; howmany = 1, which = :SR, kryldim)
    if length(model) != 0
        return eigsolve(model, howmany, which, ishermitian = true, krylovdim = kryldim)
    else
        return (missing, missing, missing)
    end
end

end
