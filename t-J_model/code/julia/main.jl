include("declarations.jl")

# Main Function
function main()
    systemSize::Int64           = 12
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 1.0
    magnonInteraction::Float64  = 0.5
    isRemovedSpinUp::Bool       = true

    if couplingJ <= 0
        println("Assumed antiferromagnetic J > 0!")
        return nothing
    end

    # calculate Heisenberg ground state energy, vector and index of magnetization subspace
    groundState::HeisenbergState = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction) #, true)

    # construct momentum space
    momentumBasis::Basis = getMomentumSpace(systemSize, groundState.magnetizationIndex, isRemovedSpinUp)

    # annihilate electron and write result in momentum basis
    holeState::Vector{ComplexF64} = getSingleHoleState(systemSize, groundState, π/2, π/2, momentumBasis, isRemovedSpinUp)

    #
    digonalizeMomentumSubspace(systemSize, tunneling, couplingJ, magnonInteraction, 0, momentumBasis)
end

@time test = main()


# TODO:
# (1) function to get possible momentum in the system
# (dependent on the size of the system)
#
# (2) asynchronous diagonalization of momentum subspaces
# (with future possibility to calculate it in parallel)
# here it is also important to merge spectral function
# inforamtion, best approach that saves memory is to
# calculate the spectral function for every momentum
# of the hole for each subspace separatelly
# in the end one should sum up the results
# for different subspace momenta
#
# (3) this time probably I would like to plot the spectral
# function in mathematica, it is also better anyway since
# I will probably run the calaculations on ICM machine,
# Therefore I shall write the spectral function to the file,
# but this also can be done in multiple ways,
# it might be better to save just overlaps end their corresponding
# eigen energies to the file and calculate the proper sum
# in mathematica, this way I will have the external control
# over broadening and will not have to rerun the calculations
# for such ridiculous reason like getting different peak width,
# this way I will store full Green's function of a single hole
# effectively
