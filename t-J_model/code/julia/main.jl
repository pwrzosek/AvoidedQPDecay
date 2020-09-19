include("declarations.jl")

# Main Function
function main()
    systemSize::Int64 = 4
    tunneling::Float64 = 1.0
    couplingJ::Float64 = tunneling * 1.0
    magnonInteraction::Float64 = 1.0

    if couplingJ <= 0
        println("The code assumes antiferromagnetic J > 0!")
        return nothing
    end

    groundState::HeisenbergState = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction, true)

    reachableSubspace::Basis = getReachableSubspace(systemSize, groundState.magnetizationIndex)

    # holeInjectedState::Vector{Float64} = getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)
end

@time test = main()
