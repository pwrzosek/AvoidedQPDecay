function initialize()
    printWelcome()
    printGeneralEvaluation()

    global g_coordinationNumber = size(g_lattice)[1]
    global g_secondCoordinationNumber = size(g_secondNeighbours)[1]
    global g_numberOfSites = size(g_lattice)[2]
    global g_numberOfMagneticConfigurations = 2^(g_numberOfSites-1)
    global g_numberOfStates = g_numberOfSites * g_numberOfMagneticConfigurations
    global g_initialStateIndex = 0
    global g_greensFunctionTemplate = Array{Any, 1}(undef, g_numberOfSites)

    printSystemInfo()
    nothing
end

function calculateAchievableSubspace()
    printAchievableSubspaceSearchBegin()

    global (g_achievableStates, g_adjacencyList) = getAchievableStatesAndAdjacencyList(g_initialStateIndex)
    global g_numberOfAchievableStates = length(g_achievableStates)
    global g_achievableSelfAdjacency = getSelfAdjacency.(g_achievableStates) .+ getImpurityPotential.(g_achievableStates)

    printAchievableSubspaceSearchEnd()
    nothing
end

function calculateLanczos(rightStatePosition, leftStatePosition)
    printLanczosEvaluationBegin()

    diagonal = Float64[]
    offDiagonal = Float64[]
    krylovSpaceDimensions = 1

    state1 = zeros(Float64, g_numberOfAchievableStates)
    state2 = zeros(Float64, g_numberOfAchievableStates)
    if rightStatePosition == leftStatePosition
        state1[rightStatePosition] = 1.
    else
        state1[rightStatePosition] = 1. / sqrt(2.)
        state1[leftStatePosition] = 1. / sqrt(2.)
    end

    function addSpan()
        # terms proportional to J
        state2 += state1 .* g_achievableSelfAdjacency
        # terms proportional to t
        for neighbour in 1 : g_coordinationNumber::Int64
            state2[g_adjacencyList[neighbour, :]] .-= g_t .* state1
        end
        # terms proportional to t'
        for secondNeighbour in 1 : g_secondCoordinationNumber::Int64
            state2[g_adjacencyList[secondNeighbour + g_coordinationNumber::Int64, :]] .-= g_tPrime .* state1
        end
        nothing
    end

    addSpan()
    push!(diagonal, BLAS.dot(state1, state2))
    axpy!(-diagonal[end], state1, state2)
    push!(offDiagonal, BLAS.nrm2(state2))

    maxDimension = min(g_maxKrylovSpaceDimension, g_numberOfAchievableStates)

    printProgressBar(krylovSpaceDimensions, maxDimension)
    while krylovSpaceDimensions < maxDimension

        state1 .*= -offDiagonal[end]
        state2 ./= offDiagonal[end]
        state1, state2 = state2, state1

        addSpan()
        push!(diagonal, BLAS.dot(state1, state2))
        axpy!(-diagonal[end], state1, state2)
        push!(offDiagonal, BLAS.nrm2(state2))

        krylovSpaceDimensions += 1
        printProgressBar(krylovSpaceDimensions, maxDimension)
    end

    printLanczosEvaluationEnd(krylovSpaceDimensions)
    (diagonal, offDiagonal, krylovSpaceDimensions)
end

function getGreensFunctionTemplate(isRead = true, fileID = string(g_numberOfSites, "_L=", g_maxKrylovSpaceDimension, "_J=", g_j, "_t'=", g_tPrime), directory = "gftmp")
    path = string(pwd(), "/", directory, "/")
    filePole = string(path, fileID, "_pole.txt")
    fileResidue = string(path, fileID, "_residue.txt")
    if isfile(filePole) && isfile(fileResidue) && isRead
        if Threads.threadid() == 1
            println("  Reading Data Files ::")
            println("    ", fileID, "_*.txt", "\n")
        end
        readLehmansRepresentation(fileID, directory)
    else
        calculateAchievableSubspace()
        calculateGreensFunctionTemplate()
        writeLehmansRepresentation(fileID, directory)
    end
end

function calculateGreensFunctionTemplate(isParallel = true)
    translationIndices = findall(cfg -> cfg == 0, last.(getStateConfiguration.(g_achievableStates)))
    translationStates = g_achievableStates[translationIndices]
    nTranslations = length(translationIndices)
    function lanczosDiagonalizeRealSpace(jobIndex)
        index = translationIndices[jobIndex]
        state = translationStates[jobIndex]
        holePosition = getHolePosition(state)
        diagonal, offDiagonal, size = calculateLanczos(g_initialStateIndex + 1, index)
        tridiagonalMatrix = SymTridiagonal(diagonal, offDiagonal[1:(end-1)])
        eigenvalues, eigenvectors = eigen(tridiagonalMatrix::SymTridiagonal{Float64,Array{Float64,1}})
        pole = eigenvalues
        residue = abs.(eigenvectors[1, :]).^2
        g_greensFunctionTemplate[holePosition + 1] = LehmansRepresentation(pole, residue)
        nothing
    end
    # calculation in parallel or in a loop based on arg g_isParallel
    if g_isParallel
        Threads.@threads for jobIndex in 1 : nTranslations
            lanczosDiagonalizeRealSpace(jobIndex)
        end
    else
        for jobIndex in 1 : nTranslations
            lanczosDiagonalizeRealSpace(jobIndex)
        end
    end
    # not every translation is possible -> catches impossible translations and adds missing LehmansRepresentation structures
    for site in 1 : g_numberOfSites
        try
            g_greensFunctionTemplate[site]
        catch
            g_greensFunctionTemplate[site] = LehmansRepresentation([0.0], [0.0])
        end
    end
    # concatenates vertically pseudo-local greens function calculations
    # more details in: Many-Body Methods for Real Materials (Pavarini, Koch, Zhang) @ page 7.25
    achievableSites = getHolePosition.(translationStates) .+ 1
    for site in achievableSites
        if site != (getHolePosition(g_initialStateIndex) + 1)
            pole = g_greensFunctionTemplate[site].m_pole
            residue = g_greensFunctionTemplate[site].m_residue
            pole = vcat(pole, g_greensFunctionTemplate[1].m_pole)
            residue = vcat(residue, -g_greensFunctionTemplate[1].m_residue)
            g_greensFunctionTemplate[site] = LehmansRepresentation(pole, residue)
        end
    end
end

function finitialize()
    printProgramEvaluationEnd()
    nothing
end

function getHolePosition(stateIndex)
    div(stateIndex, g_numberOfMagneticConfigurations)
end

function getMagneticConfiguration(stateIndex, holePosition)
    magneticIndex = stateIndex - holePosition * g_numberOfMagneticConfigurations
    protected = 2^holePosition - 1
    (magneticIndex & protected) + 2 * (magneticIndex & (~protected))
end

function getStateConfiguration(stateIndex)
    holePosition = getHolePosition(stateIndex)
    (2^holePosition, getMagneticConfiguration(stateIndex, holePosition))
end

function isParticle(type, site, stateConfiguration) # allowed types: [type == 1] -> hole, [type == 2] -> magnon
    (2^site & stateConfiguration[type]) != 0
end

function getStateIndex(holePosition, stateConfiguration)
    protected = stateConfiguration[1] - 1
    magneticIndex = (stateConfiguration[2] & protected) + (stateConfiguration[2] & (~protected)) / 2
    Int64(magneticIndex + holePosition * g_numberOfMagneticConfigurations)
end

function getAdjacentStates(stateIndex)
    result = zeros(Int64, g_coordinationNumber + g_secondCoordinationNumber)
    holePosition = getHolePosition(stateIndex)::Int64
    stateConfiguration = getStateConfiguration(stateIndex)::Tuple{Int64,Int64}
    # Nearest Neighbours
    for bond in 1 : g_coordinationNumber::Int64
        newHolePosition = g_lattice[bond, holePosition + 1]
        holeConfiguration = 2^newHolePosition
        magneticConfiguration = stateConfiguration[2]
        if isParticle(2, newHolePosition, stateConfiguration) # if magnon at new hole position
            magneticConfiguration -= 2^newHolePosition        # annihilate the magnon at new hole position
        else                                                  # else
            magneticConfiguration += 2^holePosition           # create magnon at previous hole position
        end
        newStateConfiguration = (holeConfiguration, magneticConfiguration)
        result[bond] = getStateIndex(newHolePosition, newStateConfiguration)::Int64
    end
    # Second Nearest Neighbours
    for bond in 1 : g_secondCoordinationNumber::Int64
        newHolePosition = g_secondNeighbours[bond, holePosition + 1]
        holeConfiguration = 2^newHolePosition
        magneticConfiguration = stateConfiguration[2]
        if isParticle(2, newHolePosition, stateConfiguration)           # if magnon at new hole position
            magneticConfiguration += 2^holePosition - 2^newHolePosition # move magnon to previous hole position
        end
        newStateConfiguration = (holeConfiguration, magneticConfiguration)
        result[g_coordinationNumber + bond] = getStateIndex(newHolePosition, newStateConfiguration)::Int64
    end
    result # vector of size 8 (4 nearest neighbours + 4 second nearest neighbours)
end

function getSelfAdjacency(stateIndex)
    result = .0
    stateConfiguration = getStateConfiguration(stateIndex)::Tuple{Int64,Int64}
    for siteI in 0 : (g_numberOfSites::Int64 - 1)
        isHoleAtSiteI = isParticle(1, siteI, stateConfiguration)
        isMagnonAtSiteI = isParticle(2 ,siteI, stateConfiguration)
        for bond in 1 : g_coordinationNumber::Int64
            siteJ = g_lattice[bond, siteI + 1]
            isHoleAtSiteJ = isParticle(1, siteJ, stateConfiguration)
            isMagnonAtSiteJ = isParticle(2, siteJ, stateConfiguration)
            holeTerms = isHoleAtSiteI + isHoleAtSiteJ
            magnonTerms = isMagnonAtSiteI + isMagnonAtSiteJ
            #holeHoleTerms = isHoleAtSiteI & isHoleAtSiteJ # assumptions: single hole
            magnonMagnonTerms = isMagnonAtSiteI & isMagnonAtSiteJ
            holeMagnonTerms = (isHoleAtSiteI & isMagnonAtSiteJ) + (isHoleAtSiteJ & isMagnonAtSiteI)
            result += holeTerms + magnonTerms - 2 * magnonMagnonTerms - holeMagnonTerms #- holeHoleTerms # assumptions: single hole
        end
    end
    .25 * g_j * result # half from the Hamiltonian, second half because of a bond double counting
end

# TODO check for different potentials than point potential: is lattice numeration followed by impurity numeration? [solution: define 1D like vector for impurity, than shape will follow sites numeration]
function getImpurityPotential(stateIndex)
    g_impurityScale * g_impurityShape[getHolePosition(stateIndex) + 1]
end

function getAchievableStatesAndAdjacencyList(initialStateIndex)
    achievableStates = Int64[initialStateIndex]
    adjacencyList = Int64[]
    isAchievable = zeros(Bool, g_numberOfStates)
    isAchievable[initialStateIndex + 1] = true
    for state in achievableStates
        for adjacentState in getAdjacentStates(state)::Array{Int64,1}
            push!(adjacencyList, adjacentState)
            if !isAchievable[adjacentState + 1]
                isAchievable[adjacentState + 1] = true
                push!(achievableStates, adjacentState)
            end
        end
    end
    adjacencyList = reshape(adjacencyList, (g_coordinationNumber::Int64 + g_secondCoordinationNumber::Int64, div(length(adjacencyList), g_coordinationNumber::Int64 + g_secondCoordinationNumber::Int64)))
    ascendingOrder = sortperm(achievableStates)
    achievableStates = achievableStates[ascendingOrder]
    adjacencyList = adjacencyList[:, ascendingOrder]
    for (index, state) in enumerate(adjacencyList)
        adjacencyList[index] = minimum(searchsorted(achievableStates, state))
    end
    (achievableStates, adjacencyList)
end

function getGreensFunctionValue(k, w, E0 = 0.)
    result = 0.;
    for site = 1 : g_numberOfSites
        r = g_bravais[site]
        result += exp(1im * dot(k, r)) * sum(g_greensFunctionTemplate[site].m_residue ./ (E0 .+ w .+ g_iDelta .- g_greensFunctionTemplate[site].m_pole)) #/ g_numberOfSites # note: G(k,w) = sum_r exp(ikr) G(r,w) so it should not be divided by N nor sqrt(N).
    end
    result
end

function getSpectralFunctionValue(k, w, E0 = 0.)
    -imag(getGreensFunctionValue(k, w, E0)) / pi
end

function calculateLocalSpectralFunctionTemplate()
    calculateAchievableSubspace()
    diagonal, offDiagonal, krylovDimensions = calculateLanczos(g_initialStateIndex + 1, g_initialStateIndex + 1)
    tridiagonalMatrix = diagm(0 => diagonal, 1 => offDiagonal[1:(end-1)], -1 => offDiagonal[1:(end-1)])
    eigenvalues, eigenvectors = eigen(tridiagonalMatrix)
    global g_pole = eigenvalues
    global g_residue = abs.(eigenvectors[1, :]).^2
    nothing
end

function getLocalSpectralFunctionValue(w, E0 = 0.) #minimum(g_pole))
    -imag(sum(g_residue ./ (E0 .+ w .- g_pole .+ g_iDelta))) / pi
end

function getGroundStateResidue()
    gsEnergy = minimum([g_greensFunctionTemplate[it].m_pole[1] for it in 1:g_numberOfSites])
    spectralWeights = [imag(g_iDelta) * pi * getSpectralFunctionValue(k, gsEnergy) for k in g_brillouin]
    residue = maximum(spectralWeights)
    err = 10^-10
    inds = findall(x -> abs(x - residue) < err, spectralWeights)
    result = [(spectralWeights[inds[it]], g_brillouin[inds[it]]) for it in 1:length(inds)]
end

function getGlobalVariablesNames()
    result = []
    variablesList = string.(names(Main))
    for variableName in variablesList
        if length(variableName) > 1
            if variableName[1:2] == "g_"
                push!(result, variableName)
            end
        end
    end
    sort(result)
end

function printGlobalVariables()
    println.(getGlobalVariablesNames());
    nothing
end

function printProgressBar(currentIteration, maxIteration)
    if g_showProgressBar && g_showMessages && Threads.threadid() == 1
        progress = div(100 * currentIteration, maxIteration)
        dot = progress % 4
        progress = div(progress - dot, 4)
        bar = ""
        for colon in 1 : progress
            bar *= ":"
        end
        if dot > 1
            bar *= "."
        else
            bar *= " "
        end
        for space in 2 : (25 - progress)
            bar *= " "
        end
        if currentIteration < maxIteration
            print("\r", "    [", bar, "] ", round(Int, 100 * currentIteration / maxIteration), "% ")
        else
            println("\r", "    [*************************] 100%")
        end
    end
    nothing
end

function theTime()
    string(Time(now()), " (", today(), ")")
end

function printWelcome()
    println("\n", "------ Julia for Lanczos Method ------")
    println("Today: ", today(), "   Time: ", Time(now()), "\n")
    nothing
end

function printSystemInfo()
    if Threads.threadid() == 1
        println("  Parameters ::")
        println("    Hamiltonian :  t = ", g_t, ",  t' = ", g_tPrime, ",  J = ", g_j)
        println("    Lanczos Iterations : ", g_maxKrylovSpaceDimension, "\n")
        println("  System ::")
        println("    Lattice Type : ", g_latticeTypeName)
        println("    Lattice Sites : ", g_numberOfSites)
        println("    Total Number of States : ", g_numberOfStates, " (", g_numberOfSites, " x ", g_numberOfMagneticConfigurations,")\n")
        println("  Initial State ::")
        println("    Index : ", g_initialStateIndex)
        println("    Configuration : ", getStateConfiguration(g_initialStateIndex), "\n")
    end
    nothing
end

function printGeneralEvaluation()
    if !g_showMessages
        println("  IN PROCESS * ", theTime())
    end
    nothing
end

function printAchievableSubspaceSearchBegin()
    if Threads.threadid() == 1
        println("  Calculating Achievable Subspace ::")
        println("    IN PROCESS * ", theTime())
    end
    nothing
end

function printAchievableSubspaceSearchEnd()
    if Threads.threadid() == 1
        println("    FINISHED @ ", theTime())
        println("      Number of Achievable States : ", g_numberOfAchievableStates, "\n")
    end
    nothing
end

function printLanczosEvaluationBegin()
    if g_showMessages && Threads.threadid() == 1
        println("  Evaluating Lanczos Procedure ::")
        if !g_showProgressBar
            println("    IN PROCESS * ", theTime())
        end
    end
end

function printLanczosEvaluationEnd(krylovSpaceDimensions)
    if g_showMessages && Threads.threadid() == 1
        println("    FINISHED @ ", theTime(), "            ")
        println("      Krylov Space Dimensions : ", krylovSpaceDimensions, "\n")
    end
end

function printProgramEvaluationEnd()
    println("-- DONE @ ", theTime(), " --")
    nothing
end

function printStates(fileName = "")
    if length(fileName) > 0
        file = open(fileName, "w")
        print(file, "State Index S | ")
        print(file, " HJ @ S  | ")
        println(file, "Configuration of S")
        for stateIndex in 0 : (g_numberOfStates::Int64 - 1)
            print(file, @sprintf "%13.0f | " stateIndex)
            print(file, @sprintf "%6.1f   | " getSelfAdjacency(stateIndex))
            println(file, getStateConfiguration(stateIndex))
        end
        close(file)
    else
        print("State Index S | ")
        print(" HJ @ S  | ")
        println("Configuration of S")
        for stateIndex in 0 : (g_numberOfStates::Int64 - 1)
            print(@sprintf "%13.0f | " stateIndex)
            print(@sprintf "%6.1f   | " getSelfAdjacency(stateIndex))
            println(getStateConfiguration(stateIndex))
        end
    end
    nothing
end

# function getSemianalyticalResidue()
#     gsEnergy = minimum([g_greensFunctionTemplate[it].m_pole[1] for it in 1:g_numberOfSites])
#     err = imag(g_iDelta)
#     realResidue = Vector{Vector}(undef, g_numberOfSites)
#     for site in eachindex(realResidue)
#         indices = findall(pole -> abs(pole - gsEnergy) < err, g_greensFunctionTemplate[site].m_pole[:])
#         realResidue[site] = g_greensFunctionTemplate[site].m_residue[indices]
#     end
#     realResidue = sum.(realResidue)
#     momentumResidue = similar(g_brillouin, Complex{Float64})
#     for (kSite, k) in enumerate(g_brillouin)
#         for (rSite, r) in enumerate(g_bravais)
#             momentumResidue[kSite] += realResidue[rSite] * exp(1im * dot(k, r))
#         end
#     end
#     result = 0.
#     for residue in momentumResidue
#         if abs(imag(residue)) < 10^-10
#             result = max(result, real(residue))
#         end
#     end
#     result
# end

function writeArray(A, fileName = "outp.txt")
    file = open(fileName, "w")
    nRows = length(A[:,1])
    for raw in 1:nRows
        println(file, A[raw,:])
    end
    close(file)
end

function writeLocalSpectralFunction(fileName = "outp.txt", wStep = 0.01, wMin = -3, wMax = 7)
    W = [w for w in wMin:wStep:wMax]
    A = getLocalSpectralFunctionValue.(W)
    file = open(fileName, "w")
    for value in A
        println(file, value)
    end
    close(file)
end

function writeLehmansRepresentation(fileID = string(g_numberOfSites, "_L=", g_maxKrylovSpaceDimension, "_J=", g_j, "_t'=", g_tPrime), directory = "gftmp")
    if !isdir(directory)
        mkdir(directory)
    end
    path = string(pwd(), "/", directory, "/")
    poles = Vector{Vector}(undef, g_numberOfSites)
    residues = Vector{Vector}(undef, g_numberOfSites)
    for site in eachindex(g_greensFunctionTemplate)
        poles[site] = g_greensFunctionTemplate[site].m_pole
        residues[site] = g_greensFunctionTemplate[site].m_residue
    end
    fileName = string(path, fileID, "_pole.txt")
    file = open(fileName, "w")
    writedlm(file, poles, " ")
    close(file)
    fileName = string(path, fileID, "_residue.txt")
    file = open(fileName, "w")
    writedlm(file, residues, " ")
    close(file)
end

function readLehmansRepresentation(fileID = string(g_numberOfSites, "_L=", g_maxKrylovSpaceDimension, "_J=", g_j, "_t'=", g_tPrime), directory = "gftmp")
    path = string(pwd(), "/", directory, "/")
    fileName = string(path, fileID, "_pole.txt")
    file = open(fileName, "r")
    poles = readdlm(file)
    close(file)
    fileName = string(path, fileID, "_residue.txt")
    file = open(fileName, "r")
    residues = readdlm(file)
    close(file)
    for site in eachindex(g_greensFunctionTemplate)
        g_greensFunctionTemplate[site] = LehmansRepresentation(convert.(Float64, poles[site, findall(x -> typeof(x) != SubString{String}, poles[site, :])]), convert.(Float64, residues[site, findall(x -> typeof(x) != SubString{String}, residues[site, :])]))
    end
    nothing
end

function drawLocalSpectralFunction(wStep = 0.01, wMin = -6, wMax = 8, yMax = 2)
    W = [w for w in wMin:wStep:wMax]
    A = getLocalSpectralFunctionValue.(W)
    display(Plots.plot(W, A, ylim = (0, yMax), framestyle = :box, label = string("Local :  t' = ", g_tPrime, " t", "  &  J = ", g_j, " t"), dpi = 600))
end

function draw(case = 0, kPath = g_kPath)
    A = (k, w) -> getSpectralFunctionValue(k, w) #, g_greensFunctionTemplate[1].m_pole[1])
    W = [w for w in -6:0.001:8]
    K = [g_brillouin[k] for k in kPath]
    L = [l for l in 0 : (length(K) - 1)]
    Y = zeros(Float64, length(W), length(K))
    for (ik, k) in enumerate(K)
        for (iw, w) in enumerate(W)
            Y[iw, ik] = A(k, w)
        end
    end
    if length(kPath) == 1
        display(Plots.plot(W, Y, framestyle = :box, label = string("t' = ", g_tPrime, " t", "  &  J = ", g_j, " t"), dpi = 600)) #, ylim = (-0.02, 0.84)))  #   @  k = , g_kNames[findall(x->x==kPath, g_kPath)[1][2]] # , ylim = (0, 1)
    else
        if case == 0
            display(Plots.heatmap(L, W, Y, colorbar = true, xticks = (L, g_kNames), xrotation = 90, bottom_margin = 64px, xlabel = "k", ylabel = "w / t", title = "A(w)", dpi = 600))
        else # fix it please!
            display(Plots.contour(L, W, Y, fill = true, levels = 1000, clim = (0, 1), colorbar = false))
        end
    end
    nothing
end

function saveFigure(id = 1, ncols = 2, padding = "0", directory = "plot_files",  addJ = true, addTPrime = true)
    if !isdir(directory)
        mkdir(directory)
    end
    path = string(pwd(), "/", directory, "/")
    idStr = string(lpad(id, ncols, padding), "_")
    filePathAndName = string(path, idStr, g_numberOfSites, "sites")
    if addJ
            filePathAndName = string(filePathAndName, "_J=", g_j)
    end
    if addTPrime
        filePathAndName = string(filePathAndName,  "_t'=", g_tPrime) #, "_k=", "_d=", imag(g_iDelta))) # "_V0=", g_impurityShape[1]*g_impurityScale
    end
    png(filePathAndName)
end
