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
    println(basis)
    # model = makeModel(system)
    # calculate()
    # saveResult()
end

"Reads `input.json` file and returns `info` structure with input data. It assumes `input.json` is located in the current working directory."
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

"""
function makeBasis(system::System)::Basis
    if system.magnetization < 1 || system.magnetization > system.size + 1
        error("Wrong magnetization sector in the input file!")
    end

    ### get number of spins up
    nSpinsUp::Int64 = system.magnetization - 1
    ### note: in the code spin up -> 1, spin down -> 0

    ### calculate basis size
    basisSize::Int64 = binomial(system.size, nSpinsUp)

    ### get first (lowest) state
    ### note: 1 << n == 2^n, but former is faster
    state::Int64 = nSpinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (nSpinsUp - 1))

    ### initialize the basis
    basis::Basis = Basis(state => 1)

    ### iterate over basis states
    for index in 2:basisSize
        ### initialize new state
        newState = state

        count = 0
        ### loop over <i,j> site pairs (bonds)
        for i in 0:(system.size - 2)
            j = i + 1

            ### get numeric value at i and j bit positions
            iValue, jValue = 1 << i, 1 << j

            ### check if there is bit 1 at site i
            if newState & iValue > 0
                ### check if there is bit 1 at site j
                if newState & jValue > 0
                    ### clear bit at site i
                    newState &= ~iValue

                    ### raise counter of cleared bits
                    count += 1
                else
                    ### swap bits at i and j
                    newState = xor(newState, iValue + jValue)

                    ### add cleared bits to lowest positions
                    newState += count == 0 ? 0 : sum(n -> 1 << n, 0:(count-1))

                    break
                end
            end
        end

        ### overwrite state with newState
        state = newState

        ### add state to basis
        push!(basis, state => index)
    end

    return basis
end

end
