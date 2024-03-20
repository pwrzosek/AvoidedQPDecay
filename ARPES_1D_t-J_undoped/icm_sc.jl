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

struct SC
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    momentum::Vector{Float64}
    correlation::Array{Float64, 2}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
end

function run_sc(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    β = parameters.β
    J = parameters.J
    t = parameters.t

    howmany = 2

    iDelta = 0.02im
    ωRange = collect(-3:0.001:7)

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

    vals, vecs, info = hFactorization
    GSE, GSV = vals[1], vecs[1]

    println("\n", info)


    ### String Correlator
    tJEigenvalues = Array{Float64, 2}(undef, howmany, hSystem.size + 1)
    correlation = Array{Float64, 2}(undef, 4*hSystem.size+1, 4*hSystem.size+1)
    for k in [div(hSystem.size, 4)]
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

        @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input, howmany = howmany)

        vals, vecs, _ = tJFactorization

        correlation = getStringCorrelator(vecs[1], tJBasis, tJSystem)
        epsilon = 10^-12
        if abs(vals[2] - vals[1]) < epsilon
            correlation += getStringCorrelator(vecs[2], tJBasis, tJSystem)
            correlation /= 2.0
        end

        println()
    end

    return SC(n, β, J, t, [2 * div(hSystem.size, 4) / hSystem.size], correlation)
end

function getStringCorrelator(GSV, tJBasis, tJSystem)::Array{Float64, 2}
    #C(d,s) = <S_i^z S_{i+d}^z>_{{1}i,{0}i+s,{1}i+d}
    result = zeros(Float64, 4*tJSystem.size+1, 4*tJSystem.size+1)
    for (state, index) in tJBasis
        for s in -2*tJSystem.size:2*tJSystem.size
            if mod(s, tJSystem.size) != 0
                i = mod(tJSystem.size - 1 - s, tJSystem.size)
                iValue = 1 << i
                iBit = div(state & iValue, iValue)
                Szi = 0.5 - iBit
                for d in -2*tJSystem.size:2*tJSystem.size
                    if mod(d, tJSystem.size) != 0
                        j = mod(i + d, tJSystem.size)
                        jValue = 1 << j
                        jBit = div(state & jValue, jValue)
                        Szj = 0.5 - jBit
                        if s != d
                            result[d + 2*tJSystem.size + 1, s + 2*tJSystem.size + 1] += abs(GSV[index])^2 * Szi * Szj #* (-1)^d
                        end
                    end
                end
            end
        end
    end
    return result
end

function saveData(data::Vector{SC})
    scData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Float64}, Array{Float64, 2}}}}(undef, length(data))

    for it in 1:length(data)
        scData[it] = OrderedDict(
            "size" => data[it].size,
            "interaction" => data[it].interaction,
            "coupling" => data[it].coupling,
            "hopping" => data[it].hopping,
            "momentum" => data[it].momentum,
            "correlation" => data[it].correlation
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/sc", tail, ".json"), "w")
    JSON.print(file, scData, 1)
    close(file)

    return nothing
end

### System Parameters
t = 1.0
J = 0.4
nRange = [n for n in 16:2:16]
βRange = [1.0]

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

data = Vector{SC}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_sc(parameters[it])
end

saveData(data)
