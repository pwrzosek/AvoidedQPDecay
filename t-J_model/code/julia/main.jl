# Dependencies
using LinearAlgebra

include("structures.jl")
include("functions.jl")

# Main Function
function main()
    systemSize::Int64 = 2
    tunneling::Float64 = 1.0
    couplingJ::Float64 = tunneling * 0.4

    if couplingJ <= 0
        println("The code assumes antiferromagnetic J > 0!")
        return nothing
    end

    heisenbergGroundState::StateInfo = getHeisenbergGroundState(systemSize, couplingJ)

    # Note: the model Hamiltonian (regardels of magnon-magnon interactions) conserves magnetization
    # reachableSubspace = getReachableSubspace(heisenbergGroundState.magnetizationIndex)

    # holeInjectedState::Vector{Float64} = getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)
end

main()
