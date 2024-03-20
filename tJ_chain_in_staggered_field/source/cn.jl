using OrderedCollections
using LinearAlgebra
using DelimitedFiles
using JSON

### Needed due to a bug in libopenblas64_.dll in Julia-1.5.0
if Sys.iswindows()
    LinearAlgebra.BLAS.set_num_threads(1)
end

dir = Dict(
    "single_hole" => ""
)

include(dir["single_hole"] * "single_hole.jl")

struct CN
    size::Int64
    magnetization::Int64
    momentum::Int64
    coupling::Float64
    fluctuations::Float64
    magnonmagnon::Float64
    holemagnon::Float64
    hopping::Float64
    solver::Int64
    magnons:: Vector{Int64}
    coefficients::Vector{Float64}
end

struct Parameters
    n::Int64
    m::Int64
    k::Int64
    J::Float64
    α::Float64
    λ::Float64
    λh::Float64
    t::Float64
    solver::Int64
end

function run_cn(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    m = parameters.m
    k = parameters.k

    J = parameters.J
    α = parameters.α
    λ = parameters.λ
    λh = parameters.λh

    t = parameters.t

    solver = parameters.solver


    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)

    kRange = (2 / n) .* collect(0:n) ### k range in π

    cn = Vector{Float64}(undef, 2*n - 1)
    println("Evaluating k = ", k, " in 0:", n)
    p = mod(k, n)

    input = OrderedDict(
        "system size" => n,
        "magnetization sector" => m,
        "momentum sector" => k,
        "coupling constant" => J,
        "magnon fluctuations" => α,
        "magnon interaction" => λ,
        "proximity interaction" => λh,
        "hopping constant" => t,
        "solver type" => solver
    )

    @time system, basis, model, factorization = Main.SingleHole.run(input, howmany = 2, factor = true)

    vals, vecs, _ = factorization
    cn = getCn(vecs[1], basis, system)

    ### !now we may have 2-fold degeneracy (remove spin up or down)
    epsilon = 10^-12
    ΔE = vals[2] - vals[1]
    tMc = zeros(Float64, 2*system.size - 1)
    if abs(ΔE) < epsilon
        tMc = getCn(vecs[2], basis, system)
        cn = 0.5 .* (cn .+ tMc)
    end

    println()

    return CN(
        n, m, k, J, α, λ, λh, t, solver,
        collect(-(n - 1) : (n - 1)),
        cn
    )
end

function getCn(stateVector::Vector{ComplexF64}, basis, system)
    result = zeros(Float64, 2*system.size - 1)
    for m in -(system.size - 1) : (system.size - 1)
        state = zeros(Int64, system.size)
        if m < 0
            for it in (system.size + m):(system.size - 1)
                state[system.size - it] = 1
            end
        else
            for it in 1:m
                state[system.size - it] = 1
            end
        end

        binaryState = sum([state[it] * 2^(it-1) for it in 1:length(state)])
        result[system.size + m] = abs2(stateVector[basis[binaryState]])
    end

    return result
end

function saveData(data::Vector{CN})
    cnData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Int64}, Vector{Float64}}}}(undef, length(data))

    for it in 1:length(data)
        cnData[it] = OrderedDict(
            "size" => data[it].size,
            "magnetization" => data[it].magnetization,
            "momentum" => data[it].momentum,
            "coupling" => data[it].coupling,
            "fluctuations" => data[it].fluctuations,
            "magnonmagnon" => data[it].magnonmagnon,
            "holemagnon" => data[it].holemagnon,
            "hopping" => data[it].hopping,
            "solver" => data[it].solver,
            "magnons" => data[it].magnons,
            "coefficients" => data[it].coefficients
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/cn", tail, ".json"), "w")
    JSON.print(file, cnData, 1)
    close(file)

    return nothing
end

### System Parameters
n = 12
mag = 1
mom = div(n, 4)

t = 1.0
J = 0.4

α = 1.0
λ = 1.0
λh = 1.0

solver = 0

KRange = [0.01, 0.1, 0.5] .* J

if α == 0.0
    mom = 0
end

### Init parameter list
parameters = Vector{Parameters}()
for K in KRange
    Js = J + K
    λ = J / Js
    α = λ
    push!(parameters, Parameters(n, mag, mom, Js, α, λ, λh, t, solver))
end

### Main loop
data = Vector{CN}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_cn(parameters[it])
end

saveData(data)
