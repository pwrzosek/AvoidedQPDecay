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
    holes::UInt32
    magnons::UInt32
end

# Types
Basis = Vector{State}
