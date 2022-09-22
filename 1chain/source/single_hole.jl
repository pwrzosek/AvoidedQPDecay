module SingleHole

using OrderedCollections
using LinearAlgebra
using SparseArrays
using KrylovKit
using JSON

"""
`struct System` immutable structure for system input parameters:
# Fields
*   `size::Int64`: number of lattice sites in a chain
*   `momentum::Int64`: index indicating which momentum sector should be calculated. Index `momentum` corresponds to momentum `k = 2π * momentum / size`.
*   `magnetization::Int64`: index in range `0:floor(size/2)` indicating which magnetization sector should be calculated. Magnetization is taken without sign.
*   `coupling::Float64`: value of the coupling constant J. For ferromagnet `coupling < 0`, while for antiferromagnet `coupling > 0`.
*   `interaction::FLoat64`: parameter for scaling magnon-magnon interactions. For pure Heisenberg model `interaction = 1.0`.
"""
struct System
    size::Int64
    magnetization::Int64
    momentum::Int64
    coupling::Float64
    fluctuations::Float64
    magnonmagnon::Float64
    holemagnon::Float64
    hopping::Float64
    solver::Int64
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

"`Model === Union{Matrix, SparseMatrixCSC{Complex{Float64},Int64}}`"
Model = Union{Matrix, SparseMatrixCSC{Complex{Float64},Int64}}


"""
    run()

Run t-J model diagonalization procedure.
"""
function run(input::Union{Missing, OrderedDict} = missing; howmany = 3, factor = true, timed = true)
    system::System = if input === missing
        readInput()
    else
        System(
            input["system size"],
            input["magnetization sector"],
            input["momentum sector"],
            input["coupling constant"],
            input["magnon fluctuations"],
            input["magnon interaction"],
            input["proximity interaction"],
            input["hopping constant"],
            input["solver type"]
        )
    end

    if timed
        println()
        @time basis = makeBasis(system)
        @time model = makeModel(basis, system)
    else
        basis = makeBasis(system)
        model = makeModel(basis, system)
    end
    if factor
        if timed
            @time factorization = factorize(model, system, howmany = howmany, kryldim = max(30, 3*howmany)) # put 10 for DOS
        else
            factorization = factorize(model, system, howmany = howmany, kryldim = max(30, 3*howmany))
        end
        return system, basis, model, factorization
    else
        return system, basis, model
    end
end

"Read `input.json` file and return `System` structure with input data. It requires `input.json` file to be located in the current working directory."
function readInput()::System
    path = "./source/"
    input = JSON.parsefile(path * "input.json", use_mmap = false) # use_mmap = false is a workaroud to ensure one can change JSON file without restarting Julia
    return System(
        input["system size"],
        input["magnetization sector"],
        input["momentum sector"],
        input["coupling constant"],
        input["magnon fluctuations"],
        input["magnon interaction"],
        input["proximity interaction"],
        input["hopping constant"],
        input["solver type"]
    )
end

"""
    makeBasis(system::System) -> Basis

Return `Basis === Dict{Int64, Int64}` dictionary representing basis of given magnetization and momentum sector specified by `system.magnetization` and `system.momentum` respectively. Each index in `Basis` corresponds to state value, and each value in `Basis` corresponds to position in the basis.
"""
function makeBasis(system::System)::Basis
    if isodd(system.size)
        error("Requested 'system.size' is odd! Only even sizes are supported due to transformation to magnon basis.")
    end
    if iseven(system.magnetization) || system.magnetization < 0 || system.magnetization > (system.size - 1)
        error("Wrong magnetization sector in the input file! Possible subspaces are odd integers from 1 to system.size - 1")
    end

    ### Note: in the calculated basis representative momentum states
    ###       have hole placed at site = system.size (last site of the chain)

    ### set sublattice rotation masks (we rotate sublattice without the hole in representative states)
    mask = getMask(system)

    ### number of spins pointing up
    nSpinsUp::Int64 = div(system.size + system.magnetization - 1, 2)
    ### note: in the code `spin up === 1`, `spin down === 0`

    ### calculate magnetic subspace size
    subspaceSize::Int64 = binomial(system.size - 1, nSpinsUp)

    ### get first state (i.e. with lowest index in binary basis)
    ### note: `1 << n == 2^n`, but former is faster
    state::Int64 = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))

    ### initialize the basis
    basis::Basis = Basis()

    index = 0
    ### iterate over states within given magnetization subspace
    for _ in 1:subspaceSize
        ### perform sublattice rotation
        magnonState = sublatticeRotation(state, mask, system)
        ### now 1 represents magnon, 0 stands for empty site

        ### add magnon state to basis
        push!(basis, magnonState => (index += 1))
        ### get next state
        state = getNextState(state, system)
    end
    ### to introduce momentum states we need to mix subspaces
    ### with the same abs(#spins_up - #spins_down)
    nSpinsUp = system.size - 1 - nSpinsUp
    state = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))
    for _ in 1:subspaceSize
        ### perform sublattice rotation
        magnonState = sublatticeRotation(state, mask, system)
        ### now 1 represents magnon, 0 stands for empty site

        ### add magnon state to basis
        push!(basis, magnonState => (index += 1))
        ### get next state
        state = getNextState(state, system)
    end

    return basis
end

"""
    getMask(system::System) -> Int64

Return rotation mask for the so called B sublattice. There is no hole in the B sublattice in the representative momentum states forming the basis.
"""
function getMask(system::System)::Int64
    return sum(1 << k for k in 0:2:(system.size - 2))
end

"""
    sublatticeRotation(state::Int64, mask::Int64, system::System) -> Int64

Reverse bits according to mask.
"""
function sublatticeRotation(state::Int64, mask::Int64, system::System)::Int64
    return xor(state, mask)
end

"""
    getNextState(system::System, state::Int64) -> Int64

Return `newState::Int64` next to `state::Int64` within given magnetization sector specified by `system.magnetization`.
Two states, `state` and `newState`, are next to each other withing given magentization sector if and only if there is no state that has value between `state` and `newState`.
"""
function getNextState(state::Int64, system::System)::Int64
    count = 0
    for i in 0:(system.size - 3)
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
    viewState(state::Int64, system::System)

Return matrix representing magnons (1) and ampty sites (0). Note, hole is always in site = 2*system.size in magnon momentum representative states.
"""
function viewState(state::Int64, system::System)
    return transpose(digits(state, base = 2, pad = system.size))
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


"""
    act(operator::Function, state::Int64, basis::Basis, system::System) -> LinearCombination

Apply `operator` to `state` belonging to `basis` and returns `LinearCombination  === Dict{Int64, Complex{Float64}}` representing states with their coefficients.
"""
function act(operator::Function, state::Int64, basis::Basis, system::System)::LinearCombination
    return operator(state, basis, system)
end

"""
    hamiltonian(state::Int64, basis::Basis, system::System) -> LinearCombination

Apply Hamiltonian to `state` written in magnon momentum `basis` obtained for input `system` parameters. Returns `LinearCombination  === Dict{Int64, Complex{Float64}}` representing resulting states with their coefficients.
"""
function hamiltonian(state::Int64, basis::Basis, system::System)::LinearCombination
    ### initialize result as linear combination
    result = LinearCombination(fill(state, system.size + 1), zeros(Complex{Float64}, system.size + 1))

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

        ### workout magnons first
        ### loop over bonds without the hole
        for i in 1:(system.size - 2)
            j = i + 1

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
                coefficient = 0.5 * system.coupling   # there is no exponantial since distance to representative is 0

                ### put new state in linear combination
                ### and set its corresponding coefficient
                result.state[i + 1] = newState
                result.coefficient[i + 1] = coefficient * system.fluctuations
            end

            #### diagonal coefficient
            ## comment: after rotation bits represent magnons (0 -> no magnon, 1 -> magnon present)
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            magnons = 0.5 * system.coupling * (iBit + jBit)
            magnonMagnon = -1.0 * system.coupling * iBit * jBit # note: by definition there is no magnon at the position of the hole
            result.coefficient[1] += magnons + magnonMagnon * system.magnonmagnon
        end

        ### workout sites/bonds around the hole along chain 2
        i = system.size # hole position
        for j in [1, system.size - 1]
            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)

            #### diagonal coefficient
            magnons = 0.5 * system.coupling * jBit
            hole = 0.5 * system.coupling
            holeMagnon = -0.5 * system.coupling * jBit
            result.coefficient[1] += magnons + hole + holeMagnon * system.holemagnon

            #### off-diagonal terms
            newState = state
            ## if no magnon at site of arrival
            if jBit == 0
                newState = xor(state, iValue)
                ## leave magnon behind
            else
                newState = xor(state, jValue)
                ## annihilate magnon at site of arrival
            end
            ### calculate matrix coefficient
            coefficient = 1.0   # init
            if i - j > 1    # i - j > 1 -> j === system.size + 1 (we jump forward)
                ### so we translate backwards
                ### and exp is taken with positive coeffcient exp(+)
                newState = bitmov(newState, l, false, hb = highestBit, hv = highestValue)
                coefficient *= exp(ik)
            else
                ### else hole is at system.size-1 position (we jumped backward)
                ### so we translate forward and we have exp(-)
                newState = bitmov(newState, l, true, hb = highestBit, hv = highestValue)
                coefficient *= exp(-ik)
            end

            ### create a new entry in linear combination
            ### and set its corresponding coeffcient
            j = mod(j, system.size - 1)
            result.state[system.size + j] = newState
            result.coefficient[system.size + j] = coefficient * system.hopping
        end

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
    linearCombinationLength = system.size + 1
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
    if system.solver == 1
        return Matrix(dropzeros!(sparse(I, J, V, subspaceSize, subspaceSize, +)))
    else
        return dropzeros!(sparse(I, J, V, subspaceSize, subspaceSize, +))
    end
end

"""
    factorize(model::Model, [howmany = 1, which = :SR])

Compute eigenvalues (by default with smallest real part) and their corresponding eigenvectors.
"""
function factorize(model::Model, system::System; howmany = 1, which = :SR, kryldim = 30)
    if length(model) != 0
        if system.solver == 1
            fc = eigen(model)
            return (
                fc.values,
                [fc.vectors[:,j] for j in 1:length(fc.values)],
                missing
            )
        else
            return eigsolve(model, howmany, which, ishermitian = true, krylovdim = kryldim)
        end
    else
        return (missing, missing, missing)
    end
end

end
