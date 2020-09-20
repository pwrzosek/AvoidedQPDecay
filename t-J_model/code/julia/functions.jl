# ::::::::::::::::::::::::::::::::::::::::::: #
# :::::::::::::::: Functions :::::::::::::::: #
# ::::::::::::::::::::::::::::::::::::::::::: #

# function getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction, isDebug = false)::HeisenbergState
#     # construct a full basis -- we do not know if m-m interactions
#     # affect the magnetization of the ground state but they could,
#     # nevertheless regardels of m-m interactions model conserves
#     # the magnetization so we split the basis into orthogonal
#     # subspaces corresponding to different magnetization sectors
#     function constructBasis(systemSize)::Vector{Subspace}
#         result = Vector{Subspace}(undef, Int64(systemSize + 1))
#
#         # initialize each subspace with proper size
#         for index in 1:length(result)
#             subspaceSize::Int64 = binomial(systemSize, index - 1)
#             result[index] = Subspace(index, Vector{Int64}(undef, subspaceSize))
#         end
#
#         # we want to store how many states we already assigned
#         # to each subspace to know where to put the next state
#         last = zeros(Int64, systemSize + 1)
#         numberOfStates::Int64 = 2 ^ systemSize
#
#         # the rotation is necessary in order to have a proper
#         # definition of magnons in the system,
#         # we will include the rotation of the sublattice
#         # by grouping the states with respect to staggered
#         # magnetization instead of magnetization,
#         # we choose the convention that the site at the origin
#         # (i.e. the lowest bit in binary representation) and its
#         # corresponding sublattice is not subject to the rotation
#         # while the other one sublattice is subject to the rotation
#         rotationMask::Vector{Bool} = Bool.([mod(s, 2) for s in 0:(systemSize-1)])
#
#         # for each state find its subspace
#         for state in 0:(numberOfStates - 1)
#             # we start by taking bit representation of the state
#             bitState::Vector{Bool} = digits(Bool, state, base = 2, pad = systemSize)
#
#             # we calculate staggered magnetization using rotation mask
#             # and transform it to index of corresponding subspace
#             index::Int64 = sum(bitState .⊻ rotationMask) + 1
#
#             # in the end we assign the state to its position
#             # in the basis of the proper subspace
#             result[index].basis[last[index] += 1] = state
#         end
#
#         return result
#     end
#
#     # having the basis we want to calculate the matrix for each
#     # magnetization subspace
#     function calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, subspace)::Array{Float64,2}
#         # we want to be able to apply Hamiltonian to a given state,
#         # this will produce sum of its adjacent states multiplied
#         # by the transition coefficinets (including self-transition),
#         # as the input we rather take the position of the state
#         # in the ordered basis of the corresponding subspace,
#         # as the output we want to return tuple containing
#         # positions of adjacent states (indices in the subspace basis)
#         # and transition energies (coefficients of the subspace matrix)
#         function applyHamiltonian(systemSize, couplingJ, magnonInteraction, subspace, position)::Tuple{Vector{Int64}, Vector{Float64}}
#             # we start by taking bit representation of the state
#             # corresponding to given position in the subspace basis
#             bitState::Vector{Bool} = digits(Bool, subspace.basis[position], base = 2, pad = systemSize)
#
#             # a-priori we do not know how many distinct adjacent states
#             # there exists for a chosen state
#             indices = Vector{Int64}(undef, 0)
#             coefficients = Vector{Float64}(undef, 0)
#
#             # we start by pushing the position of the initial state
#             push!(indices, position)
#             push!(coefficients, 0)
#
#             transitionCoefficient = couplingJ * 0.5
#
#             # we loop over site i in a state
#             for i in 1:systemSize
#                 # we take the neighbouring site j with periodic boundaries
#                 j = mod1(i + 1, systemSize)
#
#                 # we take care of the diagonal coefficient first
#                 coefficients[1] += couplingJ * (0.5 * (bitState[i] + bitState[j]) - 0.25 - magnonInteraction * (bitState[i] * bitState[j]))
#
#                 # for off-diagonal coefficients there are three possible scenarios
#                 # (1) two magnons, (2) two empty sites or,
#                 # (3) single magnon either in site i or site j
#                 # we care about (1) and (2) only (this comes from the transformation
#                 # we use to arrive at the bosonic model, i.e. the sublattice rotation)
#                 if (bitState[i] & bitState[j]) || ~(bitState[i] || bitState[j])
#                     # in Julia all the mutable stucts (including arrays) are
#                     # passed by reference so we have to copy bitState,
#                     # note that using simply newBitState = bitState just creates
#                     # another reference to the same address in memory and
#                     # it would  end up in changing the original bitState
#                     # which we do not want to change
#                     newBitState = copy(bitState)
#
#                     # we flip pair of spins (create/annihilate pair of magnons)
#                     newBitState[i], newBitState[j] = ~newBitState[i], ~newBitState[j]
#
#                     # we calculate new state index in binary basis
#                     newState::Int64 = sum(newBitState[it] * 2^(it-1) for it in 1:systemSize)
#
#                     # we search for the new state position in the subspace basis
#                     newPosition = searchsorted(subspace.basis, newState)[1]
#
#                     # we check if the new position is already included
#                     # and we take care of incrementing coefficients
#                     isIncluded = false
#                     for it in 1:length(indices)
#                         if indices[it] == newPosition
#                             coefficients[it] += transitionCoefficient
#                             isIncluded = true
#                             break
#                         end
#                     end
#
#                     # if it is not included we include it
#                     if !isIncluded
#                         push!(indices, newPosition)
#                         push!(coefficients, transitionCoefficient)
#                     end
#                 end
#             end
#
#             return (indices, coefficients)
#         end
#
#         # we initialize the matrix with proper size
#         # and we calculate its coeffcients
#         # here calculation for each row may be done
#         # by different threads/processes/workers [parallel]
#         subspaceSize = length(subspace.basis)
#         result = zeros(Float64, subspaceSize, subspaceSize)
#         for position in 1:subspaceSize
#             indices::Vector{Int64}, coefficients::Vector{Float64} = applyHamiltonian(systemSize, couplingJ, magnonInteraction, subspace, position)
#             result[indices[:], position] .= coefficients[:]
#         end
#
#         return result
#     end
#
#     # then we want to be able to diagonalize the matrix
#     function diagonalize(basis, isCheckDegeneracy = false)
#         # we loop over subspaces and diagonalize them,
#         # for each we extract its ground state info,
#         # each subspace can be diagonalised separately
#         # by different thread/process/worker [parallel]
#         states = Vector{HeisenbergState}(undef, systemSize + 1)
#         nSubspaces = length(basis)
#         for index in 1:nSubspaces
#             factorization = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, basis[index]))
#             # we do not expect any subspace having degenerated ground state
#             # since states are sorted in ascending order with respect
#             # to the energy we take the first one for each subspace
#             states[index] = HeisenbergState(index, factorization.values[1], factorization.vectors[:, 1])
#         end
#
#         # here we should not get degeneracy for negative J
#         # but we can check
#         if isCheckDegeneracy
#             println()
#             println("> Subspace Ground States :")
#             println(
#                 ">>> ",
#                 [state.magnetizationIndex for state in states],
#                 "\n -> ",
#                 [state.energy for state in states]
#             )
#         end
#
#         return states
#     end
#
#     # having ground states for each subspace we search
#     # for the true ground state of the system
#     function getGroundState(states)
#         result::HeisenbergState = states[1]
#         for state in states
#             if state.energy < result.energy
#                 result = state
#             end
#         end
#
#         return result
#     end
#
#     # we can also check if the ground state is not
#     # degenerated within its own subspace
#     # we define function below just for debuging purposes
#     # to make sure everything looks ok
#     function d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, heisenbergState, maxCount = 2)
#         eval = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, basis[heisenbergState.magnetizationIndex])).values
#         println()
#         println("> Degeneracy Check :")
#         println(
#             ">>> ",
#             [heisenbergState.magnetizationIndex],
#             "\n -> ",
#             eval[1:min(length(eval), maxCount)]
#         )
#     end
#
#     # construct the basis,
#     basis = constructBasis(systemSize)
#
#     # diagonalize each subspace and return the ground state info
#     groundState = getGroundState(diagonalize(basis, true))
#
#     # make sure there is no degeneracy
#     # only for debuging (takes time to run diagonalisation)
#     if isDebug
#         d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, groundState)
#     end
#
#     return groundState
# end

function getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction, isDebug = false)#::HeisenbergState
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

            # we calculate magnetization and transform it
            # to the index of corresponding subspace
            index::Int64 = sum(bitState) + 1

            # in the end we assign the state to its position
            # in the basis of the proper subspace
            result[index].basis[last[index] += 1] = state
        end

        return result
    end

    # the rotation is necessary in order to have a proper
    # definition of magnons in the system,
    # it it better to store not rotated states in the basis
    # and to apply the rotation to the state only
    # when we want to act with the Hamiltonian
    # we choose the convention that the site at the origin
    # (i.e. the lowest bit in binary representation) and its
    # corresponding sublattice is not subject to the rotation
    # while the other one sublattice is subject to the rotation
    rotationMask::Int64 = sum([[mod(s, 2) for s in 0:(systemSize-1)][it] * 2^(it-1) for it in 1:systemSize])

    # having the basis we want to calculate the matrix for each
    # magnetization subspace
    function calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, subspace)::Array{Float64,2}
        # we want to be able to apply Hamiltonian to a given state,
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
            # and we also include the rotation
            bitState::Vector{Bool} = digits(Bool, subspace.basis[position] ⊻ rotationMask, base = 2, pad = systemSize)

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

                # we take care of the diagonal coefficient first
                coefficients[1] += couplingJ * (0.5 * (bitState[i] + bitState[j]) - 0.25 - magnonInteraction * (bitState[i] * bitState[j]))

                # for off-diagonal coefficients there are three possible scenarios
                # (1) two magnons, (2) two empty sites or,
                # (3) single magnon either in site i or site j
                # we care about (1) and (2) only (this comes from the transformation
                # we use to arrive at the bosonic model, i.e. the sublattice rotation)
                if (bitState[i] & bitState[j]) || ~(bitState[i] || bitState[j])
                    # in Julia all the mutable stucts (including arrays) are
                    # passed by reference so we have to copy bitState,
                    # note that using simply newBitState = bitState just creates
                    # another reference to the same address in memory and
                    # it would end up in changing the original bitState
                    # which we do not want to change
                    newBitState = copy(bitState)

                    # we flip pair of spins (create/annihilate pair of magnons)
                    newBitState[i], newBitState[j] = ~newBitState[i], ~newBitState[j]

                    # we calculate new state index in binary basis
                    # since basis is written without rotation
                    # we rotate back the rotated sublattice
                    newState::Int64 = sum(newBitState[it] * 2^(it-1) for it in 1:systemSize) ⊻ rotationMask

                    # we search for the new state position in the subspace basis
                    newPosition = searchsorted(subspace.basis, newState)[1]

                    # we check if the new position is already included
                    # and we take care of incrementing coefficients
                    isIncluded = false
                    for it in 1:length(indices)
                        if indices[it] == newPosition
                            coefficients[it] += transitionCoefficient
                            isIncluded = true
                            break
                        end
                    end

                    # if it is not included we include it
                    if !isIncluded
                        push!(indices, newPosition)
                        push!(coefficients, transitionCoefficient)
                    end
                end
            end

            return (indices, coefficients)
        end

        # we initialize the matrix with proper size
        # and we calculate its coeffcients
        # here calculation for each row may be done
        # by different threads/processes/workers [parallel]
        subspaceSize = length(subspace.basis)
        result = zeros(Float64, subspaceSize, subspaceSize)
        for position in 1:subspaceSize
            indices::Vector{Int64}, coefficients::Vector{Float64} = applyHamiltonian(systemSize, couplingJ, magnonInteraction, subspace, position)
            result[indices[:], position] .= coefficients[:]
        end

        return result
    end

    # then we want to diagonalize the matrix
    function diagonalize(basis, isCheckDegeneracy = false)::Vector{HeisenbergState}
        # we loop over subspaces and diagonalize them,
        # for each we extract its ground state info,
        # each subspace can be diagonalised separately
        # by different thread/process/worker [parallel]
        states = Vector{HeisenbergState}(undef, systemSize + 1)
        nSubspaces = length(basis)
        for index in 1:nSubspaces
            factorization = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, basis[index]))
            # we do not expect any subspace having degenerated ground state
            # since states are sorted in ascending order with respect
            # to the energy we take the first one for each subspace
            states[index] = HeisenbergState(index, factorization.values[1], factorization.vectors[:, 1])
        end

        # here we should not get degeneracy for negative J
        # but we can check
        if isCheckDegeneracy
            println()
            println("> Subspace Ground States :")
            println(
                ">>> ",
                [state.magnetizationIndex for state in states],
                "\n -> ",
                [state.energy for state in states]
            )
        end

        return states
    end

    # having ground states for each subspace we search
    # for the true ground state of the system
    function getGroundState(states)::HeisenbergState
        result::HeisenbergState = states[1]
        for state in states
            if state.energy < result.energy
                result = state
            end
        end

        return result
    end

    # we can also check if the ground state is not
    # degenerated within its own subspace
    # we define function below just for debuging purposes
    # to make sure everything looks ok
    function d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, heisenbergState, maxCount = 2)
        eval = eigen(calculateSubspaceMatrix(systemSize, couplingJ, magnonInteraction, basis[heisenbergState.magnetizationIndex])).values
        println()
        println("> Degeneracy Check :")
        println(
            ">>> ",
            [heisenbergState.magnetizationIndex],
            "\n -> ",
            eval[1:min(length(eval), maxCount)]
        )
    end

    # construct the basis,
    basis = constructBasis(systemSize)

    # diagonalize each subspace and return the ground state info
    groundState = getGroundState(diagonalize(basis, true))

    # make sure there is no degeneracy
    # only for debuging (takes time to run diagonalisation)
    if isDebug
        d_checkDegeneracy(systemSize, couplingJ, magnonInteraction, basis, groundState)
    end

    return groundState
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

# function getReachableSubspace(systemSize, magnetizationIndex, isRemovedSpinUp = true)::Basis
#     # Heisenberg model conserves magnetization,
#     # this does not change when we add a hole,
#     # therefore we just need to know index
#     # of magnetization subspace of the Heisenberg
#     # model ground state and indicator of removed
#     # spin (i.e. up or down) to construst a basis
#     # of reachable states
#
#     # we start by setting number of spins up and down
#     # magnetizationIndex [low -> high] => magnetization [high -> low]
#     nSpinsDown::Int64 = magnetizationIndex - 1
#     nSpinsUp::Int64 = systemSize - nSpinsDown
#
#     # we need all the possible magnetic configuration
#     # for specified magnetization index and system size
#     # we exactly know how many such states there is
#     nMagneticStates = binomial(systemSize, nSpinsUp)
#
#     # we initialize magnetic configurations container
#     magneticStates = Vector{UInt32}(undef, nMagneticStates)
#
#     # the rotation is necessary in order to have a proper
#     # definition of magnons in the system,
#     # we include the rotation of the sublattice
#     # by grouping the states with respect to staggered
#     # magnetization instead of magnetization,
#     # we choose the convention that the site at the origin
#     # (i.e. the lowest bit in binary representation) and its
#     # corresponding sublattice is not subject to the rotation
#     # while the other one sublattice is subject to the rotation
#     rotationMask::Vector{Bool} = Bool.([mod(s, 2) for s in 0:(systemSize-1)])
#
#     # now we search for states that belong to specified
#     # magnetic subspace including rotation of the sublattice,
#     last::Int64 = 0
#     numberOfStates::Int64 = 2 ^ systemSize
#     for state::UInt32 in 0:(numberOfStates - 1)
#         # we start by taking bit representation of the state
#         bitState::Vector{Bool} = digits(Bool, state, base = 2, pad = systemSize)
#
#         # we calculate magnetization using rotation mask
#         # and transform it to index of corresponding subspace
#         index::Int64 = sum(bitState .⊻ rotationMask) + 1
#
#         # if the state belong to desired magnetization subspace
#         # we assign it to its position in the basis
#         if index == magnetizationIndex
#             magneticStates[last += 1] = state
#         end
#     end
#
#     # !!! from now on we assume to have a single hole
#     # injected to the ground state of the Heisenberg model
#     isRemovedSpinUp ? nSpinsUp -= 1 : nSpinsDown -= 1
#
#     # we calculate a size of the subspace of reachable states
#     function getSubspaceSize(systemSize::Int64, nSpinsUp::Int64)::Int64
#         return systemSize * binomial(systemSize - 1, nSpinsUp)
#     end
#
#     subspaceSize::Int64 = getSubspaceSize(systemSize, nSpinsUp)
#
#     # we need a way to store information
#     # about positions of holes and magnons,
#     # even if we have a single hole it is more convenient
#     # to store the configuration of hole(s) and magnons
#     # than just a position of a hole and configuration
#     # of magnons,
#     # to this end we introduce structure State and define
#     # type alias Basis for vector of States (see ./declarations.jl)
#
#     # now we can initialize the subspace basis with known size
#     reachableSubspace = Basis(undef, subspaceSize)
#
#     # we want to fill the basis with proper states,
#     # we can keep the basis in order treating the hole configuration
#     # as higher bits and magnon configuration as lower bits
#     # this is a good practice that will later allow
#     # for O(log n) search time instead of O(n)
#     # we will disregard some magnetic configurations on the way
#     # since we are going to remove either spin up or down
#     # and some configurations will simply not fit
#     # thus we shall remember the position of last assignment
#     # to reachable subspace container while looping over
#     # different hole positions and magnon configurations
#     last = 0
#
#     # it's also convenient to rewrite rotatin mask into numeric value
#     numRotation::UInt32 = sum([rotationMask[it] * 2^(it-1) for it in 1:systemSize])
#
#     # we loop over possible hole positions
#     for holePosition in 1:systemSize
#         # we take binary representation of holes as hole configuration
#         holes::UInt32 = 2^(holePosition - 1)
#
#         # we loop over magnetic configurations in desired subspace
#         # that we previously calculated
#         for magnons in magneticStates
#             # we need to know the value of the spin
#             # at the site where we want to make a hole
#             # but magnetic states were subject to rotation
#             # so we have to reverse it
#             spins::UInt32 = magnons ⊻ numRotation
#
#             # we have a following convention:
#             # spin up -> 0 and spin down -> 1
#             # while for holes
#             # hole absent -> 0, hole present -> 1
#             # so we act accordingly
#             # we check if there is a correct spin at site
#             # where we want to put a hole
#             isRemovedSpinUp ? spins = spins : spins = ~spins
#             if (holes & spins) == 0
#                 # we increment last assignment position
#                 # and assign a state to reachable subspace
#                 reachableSubspace[last += 1] = State(holes, magnons)
#             end
#         end
#     end
#
#     return reachableSubspace
# end

function getMomentumSpace(systemSize, magnetizationIndex, isRemovedSpinUp = true)::Basis
    # we start by setting number of spins up and down
    # magnetizationIndex [low -> high] => magnetization [high -> low]
    nSpinsDown::Int64 = magnetizationIndex - 1
    nSpinsUp::Int64 = systemSize - nSpinsDown

    # we need all the possible magnetic configuration
    # for specified magnetization index and system size
    # we exactly know how many such states there is
    nMagneticStates = binomial(systemSize, nSpinsUp)

    # we initialize magnetic configurations container
    magneticStates = Vector{Int64}(undef, nMagneticStates)

    # now we search for states that belong to specified
    # magnetic subspace including rotation of the sublattice,
    last::Int64 = 0
    numberOfStates::Int64 = 2 ^ systemSize
    for state in 0:(numberOfStates - 1)
        # we start by taking bit representation of the state
        bitState::Vector{Bool} = digits(Bool, state, base = 2, pad = systemSize)

        # we calculate magnetization
        # and transform it to index of corresponding subspace
        index::Int64 = sum(bitState) + 1

        # if the state belong to desired magnetization subspace
        # we assign it to its position in the basis
        if index == magnetizationIndex
            magneticStates[last += 1] = state
        end
    end

    # !!! from now on we assume to have a single hole
    # injected to the ground state of the Heisenberg model
    isRemovedSpinUp ? nSpinsUp -= 1 : nSpinsDown -= 1

    # we calculate a size of the momentum subspace
    nMomentumStates::Int64 = binomial(systemSize - 1, nSpinsUp)

    # we need a way to store information
    # about positions of holes and magnons,
    # even if we have a single hole it is more convenient
    # to store the configuration of hole(s) and magnons
    # than just a position of a hole and configuration
    # of magnons,
    # to this end we introduce structure State and define
    # type alias Basis for vector of States (see ./declarations.jl)

    # now we can initialize the subspace basis with known size
    momentumSpace = Basis(undef, nMomentumStates)

    # we want to fill the basis with proper states,
    # we will disregard some magnetic configurations on the way
    # since we are going to remove either spin up or down
    # and some configurations will simply not fit
    # thus we shall remember the position of last assignment
    last = 0

    # we loop over magnetic configurations with desited magnetization
    for spins in magneticStates
        # we check if there is a correct spin at site
        # where we want to put a hole
        if (spins & 1) == (isRemovedSpinUp ? 0 : 1)
            # we increment last assignment position
            # and assign a state to reachable subspace
            momentumSpace[last += 1] = State(1, spins)
        end
    end

    return momentumSpace
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function getHoleInjectedState(groundStateVector, reachableSubspace)

end

bitmov(x::Int, m::Int, f::Bool = true; hb::Int = 2^(m-1), pb::Int = 2^m-1) = f ? 2x - div(x, hb) * pb : div(x, 2) + rem(x, 2) * hb
