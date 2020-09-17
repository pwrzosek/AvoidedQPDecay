# Dependencies
include("structures.jl")
include("functions.jl")

# Main Function
function main()
    systemSize = 8
    tunneling = 1.0
    couplingJ = tunneling * -0.4

    if couplingJ >= 0
        println("The code assumes antiferromagnetic J < 0!")
        return nothing
    end

    heisenbergGroundState::StateInfo =
        getHeisenbergGroundState(systemSize, couplingJ)

    # Note: the model Hamiltonian (regardels of magnon-magnon interactions) conserves magnetization
    reachableSubspace =
        getReachableSubspace(heisenbergGroundState.magnetizationIndex)

    holeInjectedState::Vector{Float64} =
        getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)

end
