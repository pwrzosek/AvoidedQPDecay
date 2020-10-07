include("declarations.jl")

# Main Function
function main()
    systemSize::Int64           = 16
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 1.0
    magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    for systemSize in [8, 6]

        # calculate Lehmann's representation of the Green's function
        greensFunction::Vector{Lehmann} =
            calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)

        # save to file
        save(greensFunction, systemSize, tunneling, couplingJ, magnonInteraction)

        # print summary
        summary()
    end
end

@time main()
