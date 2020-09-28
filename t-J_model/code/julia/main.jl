include("declarations.jl")

# Main Function
function main()
    systemSize::Int64           = 4
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 1.0
    magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    # calculate Lehmann's representation of the Green's function
    greensFunction::Vector{Lehmann} =
        calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)

    # save to file
    saveLehmannRepresentation(greensFunction)
end

@time main()


# TODO:
# (1) function to get possible momenta in the system
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
