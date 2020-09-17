function getHeisenbergGroundState(systemSize, couplingJ)::StateInfo
    # construct a full basis -- we do not know if m-m interactions
    # affect the magnetization of the ground state but they can
    # nevertheless model regardels of m-m interactions conserves
    # the magnetization so we split the basis into orthogonal
    # subspaces corresponding to different magnetization sectors
    function constructBasis(systemSize)::Vector{Subspace}
        result = Vector{Subspace}(undef, Int64(systemSize + 1))
        for index in 1:length(result)
            subspaceSize::Int64 = binomial(systemSize, index - 1)
            result[index] = Subspace(index, Vector{Int64}(undef, subspaceSize))
        end
        filling = zeros(Int64, systemSize + 1)
        numberOfStates::Int64 = 2 ^ systemSize
        for state in 0:(numberOfStates - 1)
            bitRepresentation = digits(state, base = 2, pad = systemSize)
            subspaceIndex = sum(bitRepresentation) + 1
            result[subspaceIndex].basis[filling[subspaceIndex] += 1] = state
        end
        result
    end

    # having the basis we want to calculate the matrix for each
    # magnetization subspace
    function calculateSubspaceMatrix(subspace)::Array{Float64,2}
        [1 2; 3 4] # TODO: write logic here
    end

    # then we want to be able to diagonalize the matrix
    function diagonalize(basis, isCheckDegeneracy = false)
    # we loop over subspaces and diagonalize them
    # for each we extract its ground state info
        states = Vector{StateInfo}(undef, systemSize + 1)
        for (index, subspace) in enumerate(basis)
            factorization = eigen(calculateSubspaceMatrix(subspace))
            # we do not expect any subspace having degenerated ground state
            # since states are sorted in ascending order with respect
            # to the energy we take the first one for each subspace
            states[index] = StateInfo(index, factorization.values[1], factorization.vectors[:, 1])
        end
        # here we should not get degeneracy for negative J
        # but we can check
        if isCheckDegeneracy
            println([state.energy for state in states])
        end
        states
    end

    # having ground states for each subspace we search
    # for true ground state of the system
    function getGroundStateInfo(states, isCheckDegeneracy = false)
        result::StateInfo = states[1]
        for state in states
            if state.energy < result.energy
                result = state
            end
        end
        # we can also double check if the ground state is not
        # degenerated within its own subspace
        if isCheckDegeneracy
            println(eigen(calculateSubspaceMatrix(basis[result.magnetizationIndex])).values[1:2])
        end
        result
    end

    # construct the basis here
    basis = constructBasis(systemSize)
    # diagonalize each subspace
    states = diagonalize(basis, true)
    # look for the ground state info
    getGroundStateInfo(states, true)
end

function getReachableSubspace(magnetizationSubspace, isRemoveSpinUp = true)

end

function getHoleInjectedState(groundStateVector, reachableSubspace)

end
