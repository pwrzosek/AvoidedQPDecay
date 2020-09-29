# Dependencies
using LinearAlgebra
using Printf
using Dates

include("functions.jl")

# Types
Basis = Vector{Int64}

# Structures
mutable struct HeisenbergState
    magnetizationIndex::Int64
    energy::Float64
    vector::Vector{Float64}
end

mutable struct Subspace
    magnetizationIndex::Int64
    basis::Basis
end

mutable struct Lehmann
    momentum::Float64
    energy::Vector{Float64}
    weight::Vector{Float64}
end
