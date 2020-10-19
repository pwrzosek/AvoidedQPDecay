function calculatePn(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp = true)::Vector{Float64}
    # construct momentum space
    momentumBasis::Basis = getMomentumSpace(systemSize, div(systemSize + 1, 2), isRemovedSpinUp)

    # we look for the ground state
    groundStateEnergy = Inf
    groundStateVector = Vector{ComplexF64}
    groundStateMomentum = 0

    for it in 1:systemSize
        subspaceMomentum = 2π * (it - 1) / systemSize

        # diagonalize momentum subspace for eigenstates and energies
        factorization = diagonalizeMomentumSubspace(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis)

        minEnergy = real(factorization.values[1])
        if minEnergy < groundStateEnergy
            groundStateEnergy = minEnergy
            groundStateVector = deepcopy(factorization.vectors[:, 1])
        end
    end

    subspaceSize = length(momentumBasis)

    rotationMask::Int64 = sum([[mod(s, 2) for s in 0:(systemSize-1)][it] * 2^(it-1) for it in 1:systemSize])

    result = zeros(Float64, div(systemSize, 2) + 1)
    for it in 1:subspaceSize
        magnonRepresentation::Vector{Bool} = digits(Bool, momentumBasis[it] ⊻ rotationMask, base = 2, pad = systemSize)

        # check chain in forward direction
        site::Int64 = 2
        while site <= div(systemSize, 2) && magnonRepresentation[site]
            site += 1
        end
        magnonChainLength::Int64 = site - 2

        # check chain in backward direction
        site = systemSize
        while site > div(systemSize, 2) && magnonRepresentation[site]
            site -= 1
        end
        magnonChainLength = max(magnonChainLength, systemSize - site)

        # all the exponents sums out and in the end
        # we can just look at the representative
        # thus following holds true
        result[magnonChainLength + 1] += abs(groundStateVector[it])^2
    end

    return result
end

function writePn(filename, data)
    file = open(filename, "w")
    writedlm(file, data)
    close(file)
end

function mainPn()
    # systemSize::Int64           = 12
    tunneling::Float64          = 1.0
    # couplingJ::Float64          = 1.0
    # magnonInteraction::Float64  = 1.0
    isRemovedSpinUp::Bool       = true

    systemSizeRange::Vector{Int64} = [L for L in 2:2:14]
    couplingJRange::Vector{Float64} = [0.4, 1.0]
    magnonInteractionRange::Vector{Float64} = [m for m in 0.0:0.1:1.0]

    for couplingJ in couplingJRange
        println(" > J = ", couplingJ)
        for magnonInteraction in magnonInteractionRange
            println(" >> m = ", magnonInteraction)
            for systemSize in systemSizeRange
                println(" >>> L = ", systemSize)
                Pn::Vector{Float64} = calculatePn(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp)
                filename = string("../../pn/Pn_L=", systemSize, "_J=", couplingJ, "_m=", magnonInteraction, ".txt")
                writePn(filename, Pn)
            end
        end
    end

    # print summary
    summary()
end
