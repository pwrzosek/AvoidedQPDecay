# Dependencies
using LinearAlgebra

include("functions.jl")

# Structures
mutable struct HeisenbergState
    magnetizationIndex::Int64
    energy::Float64
    vector::Vector{Float64}
end

mutable struct Subspace
    magnetizationIndex::Int64
    basis::Vector{Int64}
end

mutable struct State
    hole::Int64
    spins::Int64
end

# Types
Basis = Vector{State}
