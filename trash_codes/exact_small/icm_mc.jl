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

struct MC
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    momentum::Float64
    magnons:: Vector{Int64}
    coefficients::Vector{Float64}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
end

function run_mc(parameters::Parameters)
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

    ### tJ Ground State
    k = div(hSystem.size, 4) # or 0 for tJz

    mc = Vector{Float64}(undef, 2*hSystem.size - 1)
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

    @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input)

    vals = tJFactorization.values
    vecs = tJFactorization.vectors
    mc .+= getMc(vecs[:,1], tJBasis, tJSystem)
    nn = 1

    ### !now we may have 2-fold degeneracy (remove spin up or down)
    epsilon = 10^-12
    ΔE = vals[2] - vals[1]
    tMc = zeros(Float64, 2*tJSystem.size - 1)
    if abs(ΔE) < epsilon
        tMc = getMc(vecs[:,2], tJBasis, tJSystem)
        mc .+= tMc
        nn += 1
    end
    mc ./= nn

    return MC(n, β, J, t, k, collect(-(hSystem.size - 1) : (hSystem.size - 1)), mc)
end

function getMc(stateVector::Vector{ComplexF64}, basis, system)
    result = zeros(Float64, 2*system.size - 1)
    for m in -(system.size - 1) : (system.size - 1)
        state = zeros(Int64, system.size)
        if m < 0
            for it in (system.size + m):(system.size - 1)
                state[it] = 1
            end
        else
            for it in 1:m
                state[it] = 1
            end
        end
        if mod(m, 2) == 1
            state[system.size] = 1
        end
        binaryState = sum([state[it] * 2^(it-1) for it in 1:length(state)])
        result[system.size + m] = abs2(stateVector[basis[binaryState]])
    end
    return result
end

# function numberOfMagnons(state, system, excludeHole = true)
#     return digits(state, base = 2, pad = system.size)[1:end-ifelse(excludeHole, 1, 0)]
# end
#
# function longestMagnonChain(state, system, excludeHole = true)
#     magnons = digits(state, base = 2, pad = system.size)
#     mCount = [0, 0]
#     pass = [true, true]
#     offset = ifelse(excludeHole, 1, 0)
#     for it in 1:(length(magnons)-offset)
#         pass[1] &= magnons[it] == 1
#         if pass[1]
#             mCount[1] += 1
#         end
#         pass[2] = magnons[length(magnons) + 1 - offset - it] == 1
#         if pass[2]
#             mCount[2] += 1
#         end
#     end
#     return maximum(mCount)
# end


function saveData(data::Vector{MC})
    mcData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Int64}, Vector{Float64}}}}(undef, length(data))

    for it in 1:length(data)
        mcData[it] = OrderedDict(
            "size" => data[it].size,
            "interaction" => data[it].interaction,
            "coupling" => data[it].coupling,
            "hopping" => data[it].hopping,
            "momentum" => data[it].momentum,
            "magnons" => data[it].magnons,
            "coefficients" => data[it].coefficients
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/mc", tail, ".json"), "w")
    JSON.print(file, mcData, 1)
    close(file)

    return nothing
end

### System Parameters
t = 1.0
J = 0.4
nRange = [n for n in 12:2:12]
JpRange = [1, 0.1, 0.01, 10^-10] * J
βRange = [1.0 / (1.0 + (Jp / J) / 2.0) for Jp in JpRange]

# if length(ARGS) > 0
#     t = eval(Meta.parse(ARGS[1]))
# end
# if length(ARGS) > 1
#     J = eval(Meta.parse(ARGS[2]))
# end
# if length(ARGS) > 2
#     nRange = [n for n in eval(Meta.parse(ARGS[3]))]
# end
# if length(ARGS) > 3
#     βRange = [x for x in eval(Meta.parse(ARGS[4]))]
# end

parameters = Vector{Parameters}()
for n in nRange
    for it in 1:length(βRange)
        push!(parameters, Parameters(n, βRange[it], J + JpRange[it], t))
    end
end

data = Vector{MC}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_mc(parameters[it])
end

saveData(data)
