module Dos

using LinearAlgebra
using SparseArrays
# using KrylovKit

"`Model === SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

function run(ωRange::Vector{Float64}, iDelta::Complex{Float64}, tJEigenvalues::Vector{Float64})
    dos = spectralFunction(ωRange, iDelta, tJEigenvalues)
    return dos
end

"""
    Desc.
"""
function greensFunction(ω::Complex{Float64}, tJEigenvalues::Vector{Float64})::Complex{Float64}
    return sum(1 ./ (ω .- tJEigenvalues))
end

"""
    Desc.
"""
function spectralFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, tJEigenvalues::Vector{Float64})
    result = Vector{Float64}(undef, length(ωRange))
    for (it, ω) in enumerate(ωRange)
        result[it] = -imag(greensFunction(ω + iDelta, tJEigenvalues)) / π
    end
    return result * imag(iDelta)
end

end
