# Dependencies
using LinearAlgebra

include("structures.jl")
include("functions.jl")

# Main Function
function main()
    systemSize::Int64 = 6
    tunneling::Float64 = 1.0
    couplingJ::Float64 = tunneling * 1.0
    magnonInteraction::Float64 = 1.0

    if couplingJ <= 0
        println("The code assumes antiferromagnetic J > 0!")
        return nothing
    end

    heisenbergGroundState::StateInfo = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction)

    # Note: the model Hamiltonian (regardels of magnon-magnon interactions) conserves magnetization
    # reachableSubspace = getReachableSubspace(heisenbergGroundState.magnetizationIndex)

    # holeInjectedState::Vector{Float64} = getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)

    println()
end

@time test = main()
