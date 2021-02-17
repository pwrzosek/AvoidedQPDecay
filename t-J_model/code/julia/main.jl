include("declarations.jl")

# Main Function
function main(inpSystemSize::Int64)
    systemSize::Int64           = inpSystemSize
    tunneling::Float64          = 1.0
    couplingJ::Float64          = 0.4
    magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    # calculate Lehmann's representation of the Green's function
    greensFunction::Vector{Lehmann} =
        calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)

    # save to file
    save(greensFunction, systemSize, tunneling, couplingJ, magnonInteraction)

    # print summary
    # summary()
end

@time main(4)
