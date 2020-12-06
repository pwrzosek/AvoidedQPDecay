module Heisenberg

using OrderedCollections
using KrylovKit
using JSON

"""
`struct System` immutable structure for system input parameters:
# Fields
*   `size::Int64`: number of lattice sites
*   `momentum::Int64`: index in range `1:size` indicating which momentum sector should be calculated. Range `1:size` corresponds to momenta `k` in `0:2π/size:(2π - 2π/size)`.
*   `magnetization::Int64`: index in range `1:(size + 1)` indicating which magnetization sector should be calculated. Range `1:(size + 1)` corresponds to magnetization in `-size/2:size/2`.
*   `coupling::Float64`: value of the coupling constant J. For ferromagnet `coupling < 0`, while for antiferromagnet `coupling > 0`.
*   `interaction::FLoat64`: parameters for scaling magnon-magnon interactions. For pure Heisenberg model `interaction = 1.0`.
"""
struct System
    size::Int64
    momentum::Int64
    magnetization::Int64
    coupling::Float64
    interaction::Float64
end

"`Basis === Dict{Int64, Int64}`"
Basis = OrderedDict{Int64, Int64}

"""
    run()

Runs Heisenberg model diagonalization procedure.
"""
function run()
    system::System = readInput()
    basis::Basis = makeBasis(system)
return basis
    # model = makeModel(system)
    # calculate()
    # saveResult()
end

"Reads `input.json` file and returns `System` structure with input data. It requires `input.json` file to be located in the current working directory."
function readInput()::System
    input = JSON.parsefile("input.json", use_mmap = false) # use_mmap = false is a workaroud to ensure one can change JSON file without restarting Julia
    return System(
        input["system size"],
        input["momentum sector"],
        input["magnetization sector"],
        input["coupling constant"],
        input["magnon interaction"]
    )
end

"""
    makeBasis(system::System) -> Basis

Returns `Basis === Dict{Int64, Int64}` dictionary representing basis of given magnetization and momentum sector specified by `system.magnetization` and `system.momentum` respectively. Each index in `Basis` corresponds to state value, and each value in `Basis` corresponds to position in the basis.
"""
function makeBasis(system::System)::Basis
    if system.magnetization < 1 || system.magnetization > system.size + 1
        error("Wrong magnetization sector in the input file!")
    end

    ### get number of spins up
    nSpinsUp::Int64 = system.magnetization - 1
    ### note: in the code `spin up === 1`, `spin down === 0`

    ### calculate magnetic subspace size
    subspaceSize::Int64 = binomial(system.size, nSpinsUp)

    ### get first state (i.e. with lowest index in binary basis)
    ### note: `1 << n == 2^n`, but former is faster
    state::Int64 = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))

    ### initialize the basis
    index = 0
    basis::Basis = Basis()

    ### iterate over states within given magnetization subspace
    for _ in 1:subspaceSize

        ### check if state belongs to requested momentum subspace
        if hasMomentum(system, state)
            ### pick representative state
            repState = getRepresentative(system, state)

            ### check if repState is already included
            if get(basis, repState, nothing) === nothing
                ### if not then add repState to basis
                index += 1
                push!(basis, repState => index)
            end
        end

        ### get next state
        state = getNextState(system, state)
    end

    return basis
end

"""
    hasMomentum(system::System, state::Int64) -> Bool

Returns `True` if `state` belongs to `system.momentum` subspace or returns `False` otherwise.
"""
function hasMomentum(system::System, state::Int64)::Bool
    ### initialize some parameters for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize newState
    newState = state

    ### counter for repetitions of state in translations of newState
    count = 1;
    ### note: it starts with 1 since at the begining newState == state

    ### loop over possible translations
    for _ in 2:l
        ### calculate cyclic translation of newState
        newState = bitmov(newState, l, false, hb = highestBit, hv = highestValue)

        ### raise the counter if translation goes back to state
        if newState == state
            count += 1
        end
    end

    ### state has system.momentum if system.momentum is dividable by count
    return rem(system.momentum, count) == 0
end

"""
    getRepresentative(system::System, state::Int64) -> Int64

Returns representative state `repState::Int64` of `state::Int64` within given `system::System` parameters.
`repState::Int64` is a cyclic translation of `state::Int64` that has smallest possible value.
"""
function getRepresentative(system::System, state::Int64)::Int64
    ### initialize some parameters for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize result as state
    result = state

    ### loop over translations of state
    for _ in 2:l
        ### calculate cyclic translation of state
        state = bitmov(state, l, false, hb = highestBit, hv = highestValue)

        ### replace returned state if translation has lower value
        if result > state
            result = state
        end
    end

    ### returned value is so called representative state
    ### i.e. it is the same state as initial one (with respect
    ### to cyclic translation) but it has smallest possible
    ### value within the familly of translations of initial state
    return result
end

"""
    getNextState(system::System, state::Int64) -> Int64

Returns `newState::Int64` next to `state::Int64` within given magnetization sector specified by `system.magnetization`.
Two states, `state` and `newState`, are next to each other withing given magentization sector if and only if there is no state that has value between `state` and `newState`.
"""
function getNextState(system::System, state::Int64)::Int64
    count = 0
    ### loop over <i,j> site pairs (bonds)
    for i in 0:(system.size - 2)
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
    ### if loop ends without returning the state we return -1
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

### TODO: write the model matrix / linear map

end
