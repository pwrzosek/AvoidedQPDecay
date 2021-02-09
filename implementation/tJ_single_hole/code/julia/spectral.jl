module SpectralFunction

using LinearAlgebra
using SparseArrays

"""
    Desc.
"""
mutable struct LehmansRepresentation
    pole
    residue
end

struct Krylov
    dimension
end

"`Model === SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

function run(ωRange::Vector{Float64}, initialState::Vector{Complex{Float64}}, model::Model)
    krylov = Krylov(100)
    greensFunctionTmp = getGreensFunctionTemplate(initialState, model, krylov)
end

function calculateLanczos(initialState::Vector{Complex{Float64}}, model::Model, krylov::Krylov)
    diagonal = Float64[]
    offDiagonal = Float64[]
    krylovSpaceDimensions = 1

    ### normalization of the inital state
    initialNorm = norm(initialState)

    state1 = initialState / initialNorm
    state2 = model * state1

    push!(diagonal, real(dot(state1, state2)))
    axpy!(-diagonal[end], state1, state2)
    push!(offDiagonal, Float64(norm(state2)))

    maxDimension = min(krylov.dimension, length(initialState))

    printProgressBar(krylovSpaceDimensions, maxDimension)
    while krylovSpaceDimensions < maxDimension

        state1 .*= -offDiagonal[end]
        state2 ./= offDiagonal[end]
        state1, state2 = state2, state1

        state2 += model * state1
        push!(diagonal, real(dot(state1, state2)))
        axpy!(-diagonal[end], state1, state2)
        push!(offDiagonal, Float64(norm(state2)))

        krylovSpaceDimensions += 1
        printProgressBar(krylovSpaceDimensions, maxDimension)
    end

    (diagonal, offDiagonal, krylovSpaceDimensions)
end

"""
    Desc.
"""
function getGreensFunctionTemplate(initialState::Vector{Complex{Float64}}, model::Model, krylov::Krylov)
    ### assumptions: initial state is normalized
    ### initial state is calculated for given k
    ### momentum of the basis subspace is ~system.momentum
    diagonal, offDiagonal, size = calculateLanczos(initialState, model::Model, krylov::Krylov)
    tridiagonalMatrix = SymTridiagonal(diagonal, offDiagonal[1:(end-1)])
    eigenvalues, eigenvectors = eigen(tridiagonalMatrix)
    pole = eigenvalues
    residue = abs.(eigenvectors[1, :]).^2
    greensFunctionTemplate = LehmansRepresentation(pole, residue)
    try
        greensFunctionTemplate
    catch
        greensFunctionTemplate = LehmansRepresentation([0.0], [0.0])
    end
    return greensFunctionTemplate
end

"""
    Desc.
"""
function spectralFunction(ωRange, iDelta, initialState, model::Model)

end

function printProgressBar(currentIteration, maxIteration)
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
    nothing
end

end
