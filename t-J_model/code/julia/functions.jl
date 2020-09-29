# ::::::::::::::::::::::::::::::::::::::::::: #
# :::::::::::::::: Functions :::::::::::::::: #
# ::::::::::::::::::::::::::::::::::::::::::: #

function getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction, isDebug = false)::HeisenbergState
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
        # here calculation for each column may be done
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
            # we do not expect any subspace having degenerated ground state,
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
    groundState = getGroundState(diagonalize(basis)) #, true))

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

function getMagneticSubspace(systemSize, magnetizationIndex)::Basis
    # we start by setting number of spins up and down
    # magnetizationIndex [low -> high] => magnetization [high -> low]
    nSpinsDown::Int64 = magnetizationIndex - 1
    nSpinsUp::Int64 = systemSize - nSpinsDown

    # we need all the possible magnetic configurations
    # for specified magnetization index and system size
    # we exactly know how many such states there is
    nMagneticStates = binomial(systemSize, nSpinsUp)

    # we initialize magnetic configurations container
    result = Basis(undef, nMagneticStates)

    # now we search for states that belong to specified
    # magnetic subspace including rotation of the sublattice
    last::Int64 = 0
    numberOfStates::Int64 = 2 ^ systemSize
    for state in 0:(numberOfStates - 1)
        # we start by taking bit representation of the state
        bitState::Vector{Bool} = digits(Bool, state, base = 2, pad = systemSize)

        # we calculate magnetization and get
        # the index of corresponding subspace
        index::Int64 = sum(bitState) + 1

        # if the state belong to desired magnetization subspace
        # we assign it to its position in the basis
        if index == magnetizationIndex
            result[last += 1] = state
        end
    end
    result
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function getMomentumSpace(systemSize, magnetizationIndex, isRemovedSpinUp = true)::Basis
    # construct basis of the subspace for specified magnetization index
    magneticStates = getMagneticSubspace(systemSize, magnetizationIndex)

    nSpinsDown::Int64 = magnetizationIndex - 1
    nSpinsUp::Int64 = systemSize - nSpinsDown

    # !!! from now on we assume to have a single hole
    # injected to the ground state of the Heisenberg model
    isRemovedSpinUp ? nSpinsUp -= 1 : nSpinsDown -= 1

    # we calculate the size of the momentum subspace
    nMomentumStates::Int64 = binomial(systemSize - 1, nSpinsUp)

    # now we can initialize the subspace basis with known size
    momentumSpace = Basis(undef, nMomentumStates)

    # we want to fill the basis with states,
    # we remember the position of last assignment
    last = 0

    # we loop over magnetic configurations with desired magnetization
    for spins in magneticStates
        # we check if there is a correct spin at site
        # where we want to put a hole
        if (spins & 1) == (isRemovedSpinUp ? 0 : 1)
            # we increment last assignment position
            # and assign a state to momentum space
            momentumSpace[last += 1] = spins
        end
    end

    return momentumSpace
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

# cyclic bit shift for calculationg state translations
# s - spin state configuration in binary representation
# l - number of lattice sites
# f - true -> move forward (~mult by 2); false -> backward (~div by 2)
# by ~mult and ~div we mean multiplication and division with
# preservance of periodic boundary conditions within lattice size l
# hb - highest bit [optional] (for speed up put the value of 2^(l-1))
# hv - highest value [optional] (for speed put (2^l)-1)
# if speed is not crucial one can skip hb and hv
@inline bitmov(s::Int, l::Int, f::Bool = true; hb::Int = 2^(l-1), hv::Int = 2^l-1) = f ? 2s - div(s, hb) * hv : div(s, 2) + rem(s, 2) * hb

# ******************************************* #
# ******************************************* #
# ******************************************* #

function getSingleHoleState(systemSize, state, holeMomentum, subspaceMomentum, momentumBasis, isRemovedSpinUp = true)::Vector{ComplexF64}
    # state is written in the basis that we used to obtain
    # and to diagonalize the Heisenberg model Hamiltonian,
    # the momentum basis is derived from the former one,
    # all we have to do is to go through the subspace basis
    # and for each spin state find out whether the spin
    # at site denoted as holePostion contains the spin
    # equivalen to indicator of removed spin isRemovedSpinUp,
    # note: we denote spin up as 0 while spin down as 1

    # we start by constructing basis in which state is written
    magneticStates::Basis = getMagneticSubspace(systemSize, state.magnetizationIndex)
    nMagneticStates = length(magneticStates)

    # momentum space size
    nMomentumStates = length(momentumBasis)

    # we initialize the singla hole state with proper size
    result::Vector{ComplexF64} = zeros(ComplexF64, nMomentumStates)

    # we prepare some constants for faster bit shifts
    highestBit = 2^(systemSize - 1)
    highestValue = 2^systemSize - 1

    # we loop over all magnetic states
    for it in 1:nMagneticStates
        spins = magneticStates[it]

        # we want to create the hole with given momentum
        # so we have to sum over all the possible hole
        # positions with different phase factors
        # therefore we loop over lattice sites
        for site in 1:systemSize
            # we set site 1 as the origin
            r = site - 1

            # we check if there is a proper spin at position r
            if ((spins & 2^r) == 0) == isRemovedSpinUp
                # we have to translate state to its representative
                # this just mean that we shift bits in cycle until
                # the hole is located at the origin
                representative = spins
                for _ in 1:r
                    representative = bitmov(representative, systemSize, false, hb = highestBit, hv = highestValue)
                end
                index = searchsorted(momentumBasis, representative)[1]

                # the convention is that electron annihilation with momentum k
                # is a sum over r of electron annihilations at position r
                # multiplied by a phase factors of exp(ikr) [positive sign]
                # while for momentum state T|q) = exp(iq)|q)
                # thus when we shift back T^-r |q) = exp(-iqr)|q)
                # in the end we have factor of exp(i(k-q)r)
                result[index] += exp(im * (holeMomentum - subspaceMomentum) * r) * state.vector[it] / systemSize
                # note that not all r are possible (only those where the
                # spin was up) and the value of index may be the same for
                # different spin configurations thus the above expression
                # gives contribution not only for k = q
            end
        end
    end
    result
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function diagonalizeMomentumSubspace(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis)
    # we want to write down a matrix of the t-J model
    # Hamiltonian in momentum basis for given momentum
    # and diagonalize it,
    # let us start with writing the matrix first

    function calculateSubspaceMatrix(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis)
        # we will have to include sublattice rotation
        rotationMask::Int64 = sum([[mod(s, 2) for s in 0:(systemSize-1)][it] * 2^(it-1) for it in 1:systemSize])

        # some more info for faster dealing with translations
        highestBit = 2^(systemSize - 1)
        highestValue = 2^systemSize - 1

        # we need to know how to apply the model Hamiltonian
        # to each moemntum state
        # momentum state consists of a sum over translations
        # of representative state multiplied by phase factors
        # of exp(-iqr),
        # thus we need to know how to compare each of these
        # translations to representaives forming momentum basis,
        # this can be achieved by proper tracking the position
        # of the hole in each state to be able to translate
        # states back to its own representatives while keeping
        # proper phase factor in the mean time
        function applyHamiltonian(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis, position)::Tuple{Vector{Int64}, Vector{ComplexF64}}
            # we take the position in the momentum basis
            # as the input paramter,
            # as the output we return tuple of resulting
            # states positions and their corresponding adjacency
            # we start by taking bit representation of the state
            # corresponding to given position in the subspace basis
            # and we also include the rotation as the model
            # Hamiltonian shall be formally written in terms
            # of holes and magnons
            bitState::Vector{Bool} = digits(Bool, momentumBasis[position] ⊻ rotationMask, base = 2, pad = systemSize)

            # a-priori we do not know how many distinct adjacent states
            # there exists for a chosen state
            indices = Vector{Int64}(undef, 0)
            coefficients = Vector{ComplexF64}(undef, 0)

            # we start by pushing the position of the initial state
            push!(indices, position)
            push!(coefficients, 0)

            spinFlipCoefficient = couplingJ * 0.5
            tunnelingCoefficient = -1.0 * tunneling

            # the hole is at position 0
            # when it jumps it swaps the position
            # with neighbouring spin,
            # in 1D there are two neighbours
            # first at position 1
            # second at position N-1,
            # if the hole jumps to position 1
            # we have to rotate the state backward
            # to locate the hole again at position 0
            # when we search for representatives,
            # if the hole jumps to position N-1
            # we have to rotate forward,
            # each state consists of sum over different
            # hole positions in real space
            # mulitiplied by phase factors exp(-iqr),
            # it means that overlap of the particular
            # state in real space with its corresponding
            # momentum space state is just the phase factor
            # thus when we act with a hamilatonian
            # on particular representative state
            # we know that there are also all these translations
            # with proper phase factors
            # thanks to simple form of momentum states
            # we don't have to do anything when acting with
            # H_J part since it does not change the hole
            # position and therefore does result in another
            # representative state

            # we start with the site i = 1 containg the hole
            # and take into consideration its nearest neighbourhood
            for j in [2, systemSize]
                i = 1
                # we take care of the diagonal coefficient first
                # besides the usual cost of having a magnon
                # and cost of having a hole
                # there is only hole-magnon proximity effect to take into
                # account because hole-hole interaction does not occur when
                # there is only one hole in the system
                holeCost = 0.5
                magnonCost = 0.5 * bitState[j]
                holeMagnon = -0.5 * bitState[j]
                coefficients[1] += couplingJ * (holeCost + magnonCost + holeMagnon)

                # now we want to take into account hole tunneling
                newBitState = copy(bitState)

                # if the hole jumps to unoccupied lattice site
                # it should leave the magnon behind
                # otherwise the hole should annihilate the magnon
                # in the lattice site of its arrival,
                # also we do not erase the spin at site where the hole
                # is introduced so we have to take care of it
                # even if it does not play any role it still counts
                # in proper enumeration of momentum space states
                newBitState[i], newBitState[j] = ~newBitState[j], ~newBitState[i]

                # we calculate new state index in binary basis
                # since basis is written without rotation
                # we rotate back the rotated sublattice
                newState::Int64 = sum(newBitState[it] * 2^(it-1) for it in 1:systemSize) ⊻ rotationMask

                # position of the hole after jump with respect to the origin
                r = j - i

                # we translate state back to locate hole at the origin
                # note thet once we rotated sublattice back we dont have
                # to flip spins every singla translation
                for _ in 1:r
                    newState = bitmov(newState, systemSize, false, hb = highestBit, hv = highestValue)
                end

                # each state with the hole at the origin
                # is a representative state, we just need to find
                # its index in the momentum basis
                newPosition = searchsorted(momentumBasis, newState)[1]

                # we include phase factor in coefficient
                coef = tunnelingCoefficient * exp(-im * subspaceMomentum * r)

                # we check if the new position is already included
                # and we take care of incrementing coefficients
                isIncluded = false
                for it in 1:length(indices)
                    if indices[it] == newPosition
                        coefficients[it] += coef
                        isIncluded = true
                        break
                    end
                end

                # if it is not included we include it
                if !isIncluded
                    push!(indices, newPosition)
                    push!(coefficients, coef)
                end
            end

            # we loop over remaining sites i in a state
            for i in 2:(systemSize-1)
                # we take the neighbouring site j
                # (we do not have to include periodic boundaries)
                j = i + 1

                # we take care of the diagonal coefficient first
                coefficients[1] += couplingJ * (0.5 * (bitState[i] + bitState[j]) - 0.25 - magnonInteraction * (bitState[i] * bitState[j]))

                if (bitState[i] & bitState[j]) || ~(bitState[i] || bitState[j])
                    newBitState = copy(bitState)

                    # we flip pair of spins (create/annihilate pair of magnons)
                    newBitState[i], newBitState[j] = ~newBitState[i], ~newBitState[j]

                    # we calculate new state index in binary basis
                    # since basis is written without rotation
                    # we rotate back the rotated sublattice
                    newState::Int64 = sum(newBitState[it] * 2^(it-1) for it in 1:systemSize) ⊻ rotationMask

                    # we search for the new state position in the subspace basis
                    newPosition = searchsorted(momentumBasis, newState)[1]

                    # we check if the new position is already included
                    # and we take care of incrementing coefficients
                    isIncluded = false
                    for it in 1:length(indices)
                        if indices[it] == newPosition
                            coefficients[it] += spinFlipCoefficient
                            isIncluded = true
                            break
                        end
                    end

                    # if it is not included we include it
                    if !isIncluded
                        push!(indices, newPosition)
                        push!(coefficients, spinFlipCoefficient)
                    end
                end
            end

            return (indices, coefficients)
        end

        # once we know how to act with the Hamiltonian
        # we just have to loob over momentum space states
        # and calculate the matrix
        nMomentumStates = length(momentumBasis)
        result = zeros(ComplexF64, nMomentumStates, nMomentumStates)
        for position in 1:nMomentumStates
            indices::Vector{Int64}, coefficients::Vector{ComplexF64} = applyHamiltonian(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis, position)
            result[indices[:], position] .= coefficients[:]
        end

        return result
    end

    # in the end we diagonalize the matrix
    return eigen(calculateSubspaceMatrix(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis))
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function calculateGreensFunction(systemSize, tunneling, couplingJ, magnonInteraction, isRemovedSpinUp = true)::Vector{Lehmann}
    if couplingJ <= 0
        println("Assumed antiferromagnetic J > 0!")
        return nothing
    end

    # calculate Heisenberg ground state energy, vector and index of magnetization subspace
    groundState::HeisenbergState = getHeisenbergGroundState(systemSize, couplingJ, magnonInteraction) #, true)

    # construct momentum space
    momentumBasis::Basis = getMomentumSpace(systemSize, groundState.magnetizationIndex, isRemovedSpinUp)

    # Green's functions for momentum subspaces
    subspaceResolvedResult = Array{Lehmann}(undef, systemSize + 1, systemSize)
    for it in 1:systemSize
        subspaceMomentum = 2π * (it - 1) / systemSize

        # diagonalize momentum subspace for eigenstates and energies
        factorization = diagonalizeMomentumSubspace(systemSize, tunneling, couplingJ, magnonInteraction, subspaceMomentum, momentumBasis)
        subspaceSize = length(factorization.values)

        for jt in 1:(systemSize + 1)
            holeMomentum = 2π * (jt - 1) / systemSize

            # annihilate electron and write result in momentum basis
            holeState::Vector{ComplexF64} = getSingleHoleState(systemSize, groundState, holeMomentum, subspaceMomentum, momentumBasis, isRemovedSpinUp)

            # initialize containers
            subspaceResolvedResult[jt, it] =
                Lehmann(
                    holeMomentum,
                    Vector{Float64}(undef, subspaceSize),
                    Vector{Float64}(undef, subspaceSize)
                )

            # fill in energies
            subspaceResolvedResult[jt, it].energy .= real.(factorization.values)

            # calculate overlaps
            for st in 1:subspaceSize
                subspaceResolvedResult[jt, it].weight[st] = abs(dot(factorization.vectors[:, st], holeState))^2
            end
        end
    end

    # most of weights is close to 0
    # we do not want to write them to the file
    # since they carry no information
    # thus we shall filter them out
    cutoff = 10^-20

    result = Vector{Lehmann}(undef, systemSize + 1)
    for it in 1:(systemSize + 1)
        # initialize greens function container
        result[it] =
            Lehmann(
                subspaceResolvedResult[it, 1].momentum,
                Vector{Float64}(undef, 0),
                Vector{Float64}(undef, 0)
            )
        for jt in 1:systemSize
            indices = findall(x -> x > cutoff, subspaceResolvedResult[it, jt].weight)
            append!(result[it].energy, subspaceResolvedResult[it, jt].energy[indices] .- groundState.energy)
            append!(result[it].weight, subspaceResolvedResult[it, jt].weight[indices])
        end
        ordering = sortperm(result[it].energy)
        result[it].energy = result[it].energy[ordering]
        result[it].weight = result[it].weight[ordering]
    end

    return result
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function save(greensFunction::Vector{Lehmann}, systemSize, tunneling, couplingJ, magnonInteraction)
    path = "../../data/"
    filename = string("lehmann_l=", systemSize, "_t=", tunneling, "_J=", couplingJ, "_m=", magnonInteraction,"_.txt")
    file = open(string(path, filename), "w")
    tab = "    "

    len = length(greensFunction)
    write(file, string(len, "\n\n"))
    for it in 1:len
        write(file, "{\n")
        write(file, string(tab, "{", tab, @sprintf("%3.20f", greensFunction[it].momentum), tab, "}\n"))
        for jt in 1:length(greensFunction[it].energy)
            if jt == 1
                write(file, string(tab, "{", tab))
            end
                write(file, string(@sprintf("%3.20f", greensFunction[it].energy[jt])))
            if jt != length(greensFunction[it].energy)
                write(file, tab)
            else
                write(file, tab, "}\n")
            end
        end
        for jt in 1:length(greensFunction[it].weight)
            if jt == 1
                write(file, string(tab, "{", tab))
            end
                write(file, string(@sprintf("%3.20f", greensFunction[it].weight[jt])))
            if jt != length(greensFunction[it].weight)
                write(file, tab)
            else
                write(file, tab, "}\n")
            end
        end
        write(file, "}\n\n")
    end

    close(file)
end

# ******************************************* #
# ******************************************* #
# ******************************************* #

function summary()
    println()
    println("Finished : ", Date(now()), " @ ", Time(now()))
    nothing
end
