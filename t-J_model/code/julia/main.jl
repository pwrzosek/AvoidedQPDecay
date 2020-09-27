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

    # construct momentum space
    momentumBasis::Basis = getMomentumSpace(systemSize, groundState.magnetizationIndex, isRemovedSpinUp)

    # annihilate electron and write result in momentum basis
    holeState::Vector{ComplexF64} = getSingleHoleState(systemSize, groundState, π/2, π/2, momentumBasis, isRemovedSpinUp)
end

test = main()
