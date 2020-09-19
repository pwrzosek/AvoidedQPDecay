include("declarations.jl")

# Main Function
function main()
    systemSize::Int64 = 12
    tunneling::Float64 = 1.0
    couplingJ::Float64 = tunneling * 1.0
    magnonInteraction::Float64 = 1.0

    if couplingJ <= 0
        println("The code assumes antiferromagnetic J > 0!")
        return nothing
    end

    heisenbergGroundState::StateInfo = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction)

    # Note: the model Hamiltonian (regardels of magnon-magnon interactions) conserves magnetization
    reachableSubspace::Basis = getReachableSubspace(systemSize, heisenbergGroundState.magnetizationIndex)

    # holeInjectedState::Vector{Float64} = getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)
end

@time test = main()
