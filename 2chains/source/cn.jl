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
    chain1type::Int64
    chain2type::Int64
    magnonmagnon1::Float64
    magnonmagnon2::Float64
    holemagnon1::Float64
    interchain::Float64
    interchaintype::Int64
    magnonmagnoninter::Float64
    holemagnoninter::Float64
    hopping::Float64
    periodic::Int64
    solver::Int64
    magnons:: Vector{Int64}
    coefficients::Vector{Float64}
end

struct Parameters
    n::Int64
    t::Float64
    J::Float64
    K::Float64
    α1::Int64
    α2::Int64
    αi::Int64
    λ1::Float64
    λ2::Float64
    λi::Float64
    λh1::Float64
    λhi::Float64
    pbc::Int64
    solver::Int64
end

function run_cn(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    t = parameters.t
    J = parameters.J
    K = parameters.K

    α1 = parameters.α1
    α2 = parameters.α2
    αi = parameters.αi

    λ1 = parameters.λ1
    λ2 = parameters.λ2
    λi = parameters.λi

    λh1 = parameters.λh1
    λhi = parameters.λhi

    pbc = parameters.pbc
    solver = parameters.solver


    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)

    kRange = (2 / n) .* collect(0:n) ### k range in π

    ### tJ Ground State
    k = div(n, 4)

    cn = Vector{Float64}(undef, 2*n - 1)
    println("Evaluating k = ", k, " in 0:", n)
    p = mod(k, n)

    input = OrderedDict(
        "system size" => n,
        "magnetization sector" => 1,
        "momentum sector" => p,
        "coupling constant" => J,
        "primary coupling type" => α1,
        "secondary coupling type" => α2,
        "primary magnon interaction" => λ1,
        "secondary magnon interaction" => λ2,
        "primary proximity interaction" => λh1,
        "interchain coupling" => K,
        "interchain coupling type" => αi,
        "interchain magnon interaction" => λi,
        "interchain proximity interaction" => λhi,
        "hopping constant" => t,
        "periodic along Y" => pbc,
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
        n,
        input["magnetization sector"],
        input["momentum sector"],
        J,
        α1,
        α2,
        λ1,
        λ2,
        λh1,
        K,
        αi,
        λi,
        λhi,
        t,
        pbc,
        solver,
        collect(-(n - 1) : (n - 1)),
        cn
    )
end

function getCn(stateVector::Vector{ComplexF64}, basis, system)
    result = zeros(Float64, 2*system.size - 1)
    for m in -(system.size - 1) : (system.size - 1)
        state = zeros(Int64, 2*system.size)
        if m < 0
            for it in (system.size + m):(system.size - 1)
                state[2*system.size - it] = 1
            end
        else
            for it in 1:m
                state[2*system.size - it] = 1
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
            "chain1type" => data[it].chain1type,
            "chain2type" => data[it].chain2type,
            "magnonmagnon1" => data[it].magnonmagnon1,
            "magnonmagnon2" => data[it].magnonmagnon2,
            "holemagnon1" => data[it].holemagnon1,
            "interchain" => data[it].interchain,
            "interchaintype" => data[it].interchaintype,
            "magnonmagnoninter" => data[it].magnonmagnoninter,
            "holemagnoninter" => data[it].holemagnoninter,
            "hopping" => data[it].hopping,
            "periodic" => data[it].periodic,
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
t = 1.0

J = 0.4
α1 = 1
α2 = 0

λ1 = 1.0
λ2 = 1.0

λh1 = 1.0

KRange = [0.99, 0.9, 0.5] #[2 / 99, 2 / 9, 2] .* J
αi = 0

λi = 1.0
λhi = 1.0

n = 12

pbc = 0
solver = 0

### Init parameter list
parameters = Vector{Parameters}()
for K in KRange
    J = K * 0.4
    K = 2 * J * (1.0 / K - 1)
    push!(parameters, Parameters(n, t, J, K, α1, α2, αi, λ1, λ2, λi, λh1, λhi, pbc, solver))
end

### Main loop
data = Vector{CN}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_cn(parameters[it])
end

saveData(data)
