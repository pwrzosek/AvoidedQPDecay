include("declarations.jl")

# Main Function
function main()
    systemSize::Int64           = 4
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 1.0
    magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    if couplingJ <= 0
        println("Assumed antiferromagnetic J > 0!")
        return nothing
    end

    # calculate Heisenberg ground state energy, vector and index of magnetization subspace
    groundState::HeisenbergState = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction, true)

    # # construct basis of reachable states starting from the specified magnetic subspace
    # reachableSubspace::Basis = getReachableSubspace(systemSize, groundState.magnetizationIndex, isRemovedSpinUp)
    #
    # construct momentum space
    momentumSpace::Basis = getMomentumSpace(systemSize, groundState.magnetizationIndex, isRemovedSpinUp)

    # holeInjectedState::Vector{Float64} = getHoleInjectedState(heisenbergGroundState.vector, reachableSubspace)
end

@time test = main()
