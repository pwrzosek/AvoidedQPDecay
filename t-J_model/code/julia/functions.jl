function getHeisenbergGroundState(systemSize, couplingJ)::StateInfo
    # construct a full basis -- we do not know id m-m interactions
    # affect the magnetization of the ground state but they can
    # nevertheless model regardels of m-m interactions conserves
    # the magnetization so we split the basis into orthogonal
    # subspaces corresponding to different magnetization sectors
    function constructBasis(systemSize)

    end

    # having the basis we want to calculate the matrix for each
    # magnetization subspace
    function calculateSubspaceMatrix(subspace)

    end

    # then we want to be able to diagonalize the matrix
    function diagonalize(basis, isCheckDegeneracy = false)
    # we loop over subspaces and diagonalize them
    # for each we extract its ground state info
        states = Vector{StateInfo}(undef, 2 * systemSize + 1)
        for (index, subspace) in enumerate(basis)
            factorization = eigen(subspace)
            # we do not expect any subspace having degenerated ground state
            states[index] = StateInfo(factorization.vectors[:, 1], factorization.values[1], index)
        end
        # in the end we look for the subspace with lowest energy
        # here we should not get degeneracy for negative J
        # but we check can check
        if isCheckDegeneracy
            println([state.energy for state in states])
        end
        states
    end

    # we search for true ground state of the system
    function getGroundStateInfo(states, isCheckDegeneracy = false)
        result = states[1]
        for state in states
            if state.energy < result.energy
                result = state
            end
        end
        # we can also double check if the ground state is not
        # degenerated within its own subspace
        if isCheckDegeneracy
            nextEnergy = diagonalize(basis[result.magnetizationIndex]).values[2]
            println([result.energy, nextEnergy])
        end
        result
    end

    # we construct the basis here
    basis = constructBasis(systemSize)
    # we diagonalize each subspace
    states = diagonalize(basis, true)
    # and we search for true ground state
    getGroundStateInfo(states, true)
end

function getReachableSubspace(magnetizationSubspace, isRemoveSpinUp = true)

end

function getHoleInjectedState(groundStateVector, reachableSubspace)

end
