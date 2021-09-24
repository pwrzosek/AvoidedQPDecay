using OrderedCollections
using LinearAlgebra
using DelimitedFiles
using JSON

### Needed due to a bug in libopenblas64_.dll in Julia-1.5.0
if Sys.iswindows()
    LinearAlgebra.BLAS.set_num_threads(1)
end

dir = Dict(
    "heisenberg" => "heisenberg/code/julia/",
    "tJModel" => "tJ_single_hole/code/julia/"
)

include(dir["heisenberg"] * "heisenberg.jl")
include(dir["tJModel"] * "tJmodel.jl")
include(dir["tJModel"] * "spectral.jl")

struct PN
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    momentum::Vector{Float64}
    P0n::Vector{Float64}
    P1n::Vector{Vector{Float64}}
    # energy::Vector{Float64}
    # spectrum::Array{Float64, 2}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
end

function run_pn(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    β = parameters.β
    J = parameters.J
    t = parameters.t

    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)

    kRange = (2 / n) .* collect(0:n) ### k range in π
    q = 0 # in rotating momentum basis Heisenberg GS always at q = 0

    ### Heisenberg Ground State
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => q,
        "magnetization sector" => 0,
        "coupling constant" => J,
        "magnon interaction" => β
    )

    println()
    @time hSystem, hBasis, hFactorization = Main.Heisenberg.run(input)

    _, vecs, info = hFactorization

    println("\n", info)

    P0n = getPn(vecs[1], hBasis, hSystem, false)

    ### tJ Ground State
    P1n = Vector{Vector{Float64}}(undef, hSystem.size + 1)
    for k in 0:hSystem.size
        println("Evaluating k = ", k, " in 0:", hSystem.size)
        p = mod(k + q, hSystem.size)

        input = OrderedDict(
            "system size" => n,
            "momentum sector" => p,
            "magnetization sector" => 0,
            "coupling constant" => J,
            "magnon interaction" => β,
            "hopping constant" => t
        )

        @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input, howmany = 2, factor = true)

        vals, vecs, _ = tJFactorization
        P1n[k + 1] = getPn(vecs[1], tJBasis, tJSystem, true)

        ### !now we may have 2-fold degeneracy (remove spin up or down)
        epsilon = 10^-12
        ΔE = vals[2] - vals[1]
        tPn = zeros(Float64, tJSystem.size + 1)
        if abs(ΔE) < epsilon
            tPn = getPn(vecs[2], tJBasis, tJSystem)
            P1n[k + 1] = (P1n[k + 1] + tPn) / 2
        end

        println()
    end

    return PN(n, β, J, t, kRange, P0n, P1n)
end

function getPn(stateVector::Vector{ComplexF64}, basis, system, excludeHole = true, magnonCounter::Function = numberOfMagnons)
    result = zeros(Float64, system.size + 1)
    for (state, index) in basis
        magnons = magnonCounter(state, system, excludeHole)
        result[sum(magnons) + 1] += abs(stateVector[index])^2
    end
    return result
end

function numberOfMagnons(state, system, excludeHole = true)
    return digits(state, base = 2, pad = system.size)[1:end-ifelse(excludeHole, 1, 0)]
end

function longestMagnonChain(state, system, excludeHole = true)
    magnons = digits(state, base = 2, pad = system.size)
    mCount = [0, 0]
    pass = [true, true]
    offset = ifelse(excludeHole, 1, 0)
    for it in 1:(length(magnons)-offset)
        pass[1] &= magnons[it] == 1
        if pass[1]
            mCount[1] += 1
        end
        pass[2] = magnons[length(magnons) + 1 - offset - it] == 1
        if pass[2]
            mCount[2] += 1
        end
    end
    return maximum(mCount)
end


function saveData(data::Vector{PN})
    pnData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Float64}, Vector{Vector{Float64}}}}}(undef, length(data))

    for it in 1:length(data)
        pnData[it] = OrderedDict(
            "size" => data[it].size,
            "interaction" => data[it].interaction,
            "coupling" => data[it].coupling,
            "hopping" => data[it].hopping,
            "momentum" => data[it].momentum,
            "P0n" => data[it].P0n,
            "P1n" => data[it].P1n
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/pn", tail, ".json"), "w")
    JSON.print(file, pnData, 1)
    close(file)

    return nothing
end

### System Parameters
t = 1.0
J = 0.4
nRange = [n for n in 16:2:16]
βRange = [0.0, 1.0]

if length(ARGS) > 0
    t = eval(Meta.parse(ARGS[1]))
end
if length(ARGS) > 1
    J = eval(Meta.parse(ARGS[2]))
end
if length(ARGS) > 2
    nRange = [n for n in eval(Meta.parse(ARGS[3]))]
end
if length(ARGS) > 3
    βRange = [x for x in eval(Meta.parse(ARGS[4]))]
end

parameters = Vector{Parameters}()
for n in nRange
    for β in βRange
        push!(parameters, Parameters(n, β, J, t))
    end
end

data = Vector{PN}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_pn(parameters[it])
end

saveData(data)
