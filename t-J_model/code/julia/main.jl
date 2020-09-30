include("declarations.jl")

# Main Function
function main()
    systemSize::Int64           = 12
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 0.4
    magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    # for magnonInteraction in [0.4, 0.3, 0.2, 0.1, 0.0]

        # calculate Lehmann's representation of the Green's function
        greensFunction::Vector{Lehmann} =
            calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)

        # save to file
        save(greensFunction, systemSize, tunneling, couplingJ, magnonInteraction)

        # print summary
        summary()
    # end
end

@time main()
