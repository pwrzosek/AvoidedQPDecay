module SpectralFunction

using LinearAlgebra
using SparseArrays
using KrylovKit

struct Krylov
    dimension::Int64
end

"`Model === SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

function run(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, model::Model, krylovDimension::Int64 = 400)
    krylov = Krylov(krylovDimension)
    diagonal, offDiagonal, size = calculateLanczos(initialState, model::Model, krylov::Krylov)
    spectrum = spectralFunction(ωRange, iDelta, initialState, diagonal, offDiagonal, size)
    return spectrum
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

    # printProgressBar(krylovSpaceDimensions, maxDimension)
    while krylovSpaceDimensions < maxDimension

        state1 .*= -offDiagonal[end]
        state2 ./= offDiagonal[end]
        state1, state2 = state2, state1

        state2 += model * state1
        push!(diagonal, real(dot(state1, state2)))
        axpy!(-diagonal[end], state1, state2)
        push!(offDiagonal, Float64(norm(state2)))

        krylovSpaceDimensions += 1
        # printProgressBar(krylovSpaceDimensions, maxDimension)
    end

    (diagonal, offDiagonal[1:(end-1)], krylovSpaceDimensions)
end


"""
    Desc.
"""
function greensFunction(ω, initialState, diagonal, offDiagonal, size)
    result = ω - diagonal[size]
    for it in (size - 1):-1:1
        result = ω - diagonal[it] - offDiagonal[it]^2 / result
    end
    return dot(initialState, initialState) / result
end

"""
    Desc.
"""
function spectralFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, diagonal::Vector{Float64}, offDiagonal::Vector{Float64}, size::Int64)
    result = Vector{Float64}(undef, length(ωRange))
    for (it, ω) in enumerate(ωRange)
        result[it] = -imag(greensFunction(ω + iDelta, initialState, diagonal, offDiagonal, size)) / π
    end
    return result
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
        println("\r", "    [:::::::::::::::::::::::::] 100%")
    end
    nothing
end

end
