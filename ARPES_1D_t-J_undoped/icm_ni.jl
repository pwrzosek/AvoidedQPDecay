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

struct NI
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    anisotropy::Float64
    coefficients::Vector{Float64}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
    α::Float64
end

function run_ni(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    β = parameters.β
    J = parameters.J
    t = parameters.t
    α = parameters.α

    q = 0 # in rotating momentum basis Heisenberg GS always at q = 0

    ### Heisenberg Ground State
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => q,
        "magnetization sector" => 0,
        "coupling constant" => J,
        "anisotropy" => α,
        "magnon interaction" => β
    )

    println()
    @time hSystem, hBasis, hFactorization = Main.Heisenberg.run(input)

    _, vecs, info = hFactorization

    println("\n", info)

    ni = Vector{Float64}(undef, hSystem.size)
    ni = getNi(vecs[1], hBasis, hSystem)

    println()

    return NI(n, β, J, t, α, ni)
end

function getNi(stateVector::Vector{ComplexF64}, basis, system)
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    function makeCase(magnons)
        size = length(magnons)
        if sum(magnons) > size / 2
            return ones(size) - magnons
        else
            return magnons
        end
        # return magnons
    end

    result = zeros(Float64, system.size)
    for (state, index) in basis
        magnons = digits(state, base = 2, pad = system.size)

        ni = abs2(stateVector[index]) * makeCase(magnons)
        periodicity::Int64 = 1
        stateTranslation = state
        while state != (stateTranslation = Main.Heisenberg.bitmov(stateTranslation, l, false, hb = highestBit, hv = highestValue))
            magnons = digits(stateTranslation, base = 2, pad = system.size)
            ni += abs2(stateVector[index]) * makeCase(magnons)
            periodicity += 1
        end

        result += ni / periodicity

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


function saveData(data::Vector{NI})
    niData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Int64}, Vector{Float64}}}}(undef, length(data))

    for it in 1:length(data)
        niData[it] = OrderedDict(
            "size" => data[it].size,
            "interaction" => data[it].interaction,
            "coupling" => data[it].coupling,
            "hopping" => data[it].hopping,
            "anisotropy" => data[it].anisotropy,
            "coefficients" => data[it].coefficients
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/ni", tail, ".json"), "w")
    JSON.print(file, niData, 1)
    close(file)

    return nothing
end

### System Parameters
t = 1.0
J = 0.4
nRange = [n for n in 16:2:16]
αRange = [α for α in 0:0.1:4]
βRange = [1.0]

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
    for β in βRange
        for α in αRange
            push!(parameters, Parameters(n, β, J, t, α))
        end
    end
end

data = Vector{NI}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_ni(parameters[it])
end

saveData(data)
