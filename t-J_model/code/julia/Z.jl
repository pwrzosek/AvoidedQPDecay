function calculateZ(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp = true)
    # calculate Lehmann's representation of the Green's function
    greensFunction::Vector{Lehmann} =
        calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)

    result::Array{Float64,2} = Array{Float64,2}(undef, systemSize + 1, 2)
    for it in 1:(systemSize + 1)
        result[it, 1] = greensFunction[it].momentum
        result[it, 2] = greensFunction[it].weight[1]
    end

    return result
end

function writeZ(filename, data)
    file = open(filename, "w")
    writedlm(file, data)
    close(file)
end

function mainZ()
    # systemSize::Int64           = 12
    tunneling::Float64          = 1.0
    # couplingJ::Float64          = 1.0
    # magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    systemSizeRange::Vector{Int64} = [L for L in 2:2:12]
    couplingJRange::Vector{Float64} = [1.0]
    magnonInteractionRange::Vector{Float64} = [m for m in 0.0:0.1:1.0]

    for couplingJ in couplingJRange
        println(" > J = ", couplingJ)
        for magnonInteraction in magnonInteractionRange
            println(" >> m = ", magnonInteraction)
            for systemSize in systemSizeRange
                println(" >>> L = ", systemSize)
                Z::Array{Float64,2} = calculateZ(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)
                filename = string("../../z/z_L=", systemSize, "_J=", couplingJ, "_m=", magnonInteraction, ".txt")
                writeZ(filename, Z)
            end
        end
    end

    # print summary
    # summary()

end
