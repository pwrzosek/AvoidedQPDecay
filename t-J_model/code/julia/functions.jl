function getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction)::StateInfo
    # construct a full basis -- we do not know if m-m interactions
    # affect the magnetization of the ground state but they could,
    # nevertheless regardels of m-m interactions model conserves
    # the magnetization so we split the basis into orthogonal
    # subspaces corresponding to different magnetization sectors
    function constructBasis(systemSize)::Vector{Subspace}
        result = Vector{Subspace}(undef, Int64(systemSize + 1))

        # initialize each subspace with proper size
        for index in 1:length(result)
            subspaceSize::Int64 = binomial(systemSize, index - 1)
            result[index] = Subspace(index, Vector{Int64}(undef, subspaceSize))
        end

        # we want to store how many states we already assigned
        # to each subspace to know where to put the next state
        last = zeros(Int64, systemSize + 1)
        numberOfStates::Int64 = 2 ^ systemSize

        # for each state find its subspace
        for state in 0:(numberOfStates - 1)
            # we start by taking bit representation of the state
            bitState::Vector{Bool} = digits(Bool, state, base = 2, pad = systemSize)

            # we will include the rotation of the sublattice
            # by grouping the states with respect to staggered
            # magnetization instead of magnetization,
            # the rotation is necessary in order to have a proper
            # definition of magnons in the system,
            # we choose the convention that the site at origin
            # (i.e. the lowest bit in binary representation) and its
            # corresponding sublattice is not subject to the rotation
            # while the other one sublattice is subject to the rotation
            rotationMask::Vector{Bool} = Bool.([mod(s, 2) for s in 0:(systemSize-1)])

            # we calculate staggered magnetization using rotation mask
            # and transform it to index of corresponding subspace
            index::Int64 = sum(bitState .⊻ rotationMask) + 1

            # in the end we assign the state to its position
            # in the basis of the proper subspace
            result[index].basis[last[index] += 1] = state
        end

        return result
    end

    # having the basis we want to calculate the matrix for each
    # magnetization subspace
    function calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, subspace)::Array{Float64,2}
        # we want to be able to apply hamiltonian to a given state,
        # this will produce sum of its adjacent states multiplied
        # by the transition coefficinets (including self-transition),
        # as the input we rather take the position of the state
        # in the ordered basis of the corresponding subspace,
        # as the output we want to return tuple containing
        # positions of adjacent states (indices in the subspace basis)
        # and transition energies (coefficients of the subspace matrix)
        function applyHamiltonian(systemSize, couplingJ, magnonInteraction, subspace, position)::Tuple{Vector{Int64}, Vector{Float64}}
            # we start by taking bit representation of the state
            # corresponding to given position in the subspace basis
            bitState::Vector{Bool} = digits(Bool, subspace.basis[position], base = 2, pad = systemSize)

            # a-priori we do not know how many distinct adjacent states
            # there exists for a chosen state
            indices = Vector{Int64}(undef, 0)
            coefficients = Vector{Float64}(undef, 0)

            # we start by pushing the position of the initial state
            push!(indices, position)
            push!(coefficients, 0)

            transitionCoefficient = couplingJ * 0.5

            # we loop over site i in a state
            for i in 1:systemSize
                # we take the neighbouring site j with periodic boundaries
                j = mod1(i + 1, systemSize)

                # we take care of the diagonal cefficient first
                coefficients[1] += couplingJ * (0.5 * (bitState[i] + bitState[j]) - 0.25 - magnonInteraction * (bitState[i] * bitState[j]))

                # for off-diagonal coefficients there are three possible scenarios
                # (1) two magnons, (2) two empty sites or,
                # (3) single magnon either in site i or site j
                # we care about (1) and (2) only (this comes from the transformation
                # we use to arrive at the bosonic model, i.e. sublattice rotation
                if (bitState[i] & bitState[j]) || ~(bitState[i] || bitState[j])
                    # in Julia all the mutable stucts (including arrays) are
                    # passed by reference so we have to copy bitState,
                    # note that using simply newBitState = bitState just creats
                    # another reference to the same address in memory and
                    # it would  end up in changing the original bitState
                    # which we do not want to change
                    newBitState = copy(bitState)

                    # we flip pair of spins (create/annihilate pair of magnons)
                    newBitState[i], newBitState[j] = ~newBitState[i], ~newBitState[j]

                    # we calculate new state index in binary basis
                    newState::Int64 = sum(newBitState[k] * 2^(k-1) for k in 1:systemSize)

                    # we search for the new state position in the subspace basis
                    newPosition = searchsorted(subspace.basis, newState)[1]

                    # we check if the new position is already included
                    # we include it if it is not and we take care of
                    # incrementing transition coefficients
                    isIncluded = false
                    for it in 1:length(indices)
                        if indices[it] == newPosition
                            coefficients[it] += transitionCoefficient
                            isIncluded = true
                            break
                        end
                    end
                    if !isIncluded
                        push!(indices, newPosition)
                        push!(coefficients, transitionCoefficient)
                    end
                end
            end

            return (indices, coefficients)
        end

        subspaceSize = length(subspace.basis)

        # we initialize the matrix with proper size
        # and we calculate its coeffcients
        # here calculation for each column may be done
        # in async way by different threads/processes/workers
        result = zeros(Float64, subspaceSize, subspaceSize)
        for position in 1:subspaceSize
            indices::Vector{Int64}, coefficients::Vector{Float64} = applyHamiltonian(systemSize, couplingJ, magnonInteraction, subspace, position)
            result[indices[:], position] .= coefficients[:]
        end

        return result
    end

    # then we want to be able to diagonalize the matrix
    function diagonalize(basis, isCheckDegeneracy = false)
        # we loop over subspaces and diagonalize them
        # for each we extract its ground state info
        states = Vector{StateInfo}(undef, systemSize + 1)
        for (index, subspace) in enumerate(basis)
            factorization = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, subspace))
            # we do not expect any subspace having degenerated ground state
            # since states are sorted in ascending order with respect
            # to the energy we take the first one for each subspace
            states[index] = StateInfo(index, factorization.values[1], factorization.vectors[:, 1])
        end

        # here we should not get degeneracy for negative J
        # but we can check
        if isCheckDegeneracy
            println()
            println("> Subspace Ground States :")
            println(
                ">>> ",
                [state.magnetizationIndex for state in states],
                " -> ",
                [state.energy for state in states]
            )
        end

        return states
    end

    # having ground states for each subspace we search
    # the true ground state of the system
    function getGroundStateInfo(states, isCheckDegeneracy = false)
        result::StateInfo = states[1]
        for state in states
            if state.energy < result.energy
                result = state
            end
        end

        return result
    end

    # we can also check if the ground state is not
    # degenerated within its own subspace
    # we define below function just for debuging purpose
    # to make sure everything looks ok
    function d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, stateInfo, maxCount = 2)
        eval = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, basis[stateInfo.magnetizationIndex])).values
        println()
        println("> Degeneracy Check :")
        println(
            ">>> ",
            [stateInfo.magnetizationIndex],
            " -> ",
            eval[1:min(length(eval), maxCount)]
        )
    end

    # construct the basis,
    basis = constructBasis(systemSize)

    # diagonalize each subspac and return the ground state info
    stateInfo = getGroundStateInfo(diagonalize(basis, true), true)

    # we make sure there is no degeneracy
    d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, stateInfo)

    return stateInfo
end

function getReachableSubspace(magnetizationSubspace, isRemoveSpinUp = true)

end

function getHoleInjectedState(groundStateVector, reachableSubspace)

end
