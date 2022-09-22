using OrderedCollections
using LinearAlgebra
using DelimitedFiles
using JSON

### Needed due to a bug in libopenblas64_.dll in Julia-1.5.0
if Sys.iswindows()
    LinearAlgebra.BLAS.set_num_threads(1)
end

dir = Dict(
    "half_filling" => "half_filling/code/julia/",
    "single_hole" => "single_hole/code/julia/"
)

include(dir["half_filling"] * "half_filling.jl")
include(dir["single_hole"] * "single_hole.jl")
# include(dir["tJModel"] * "spectral.jl")

struct MC
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    momentum::Float64
    intrachain::Float64
    isotropy::Float64
    magnons:: Vector{Int64}
    coefficients::Vector{Float64}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
    K::Float64
    α::Float64
end

function run_mc(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    β = parameters.β
    J = parameters.J
    t = parameters.t
    K = parameters.K
    α = parameters.α


    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)

    kRange = (2 / n) .* collect(0:n) ### k range in π
    q = 0 # in rotating momentum basis Heisenberg GS always at q = 0

    ### Half-Filling Ground State
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => q,
        "magnetization sector" => 0,
        "coupling constant" => J,
        "magnon interaction" => β,
        "intrachain coupling" => K,
        "intrachain isotropy" => α
    )

    println()
    @time hSystem, hBasis, hFactorization = Main.HalfFilling.run(input)

    _, vecs, info = hFactorization

    println("\n", info)

    ### tJ Ground State
    k = div(hSystem.size, 4)

    mc = Vector{Float64}(undef, 2*hSystem.size - 1)
    println("Evaluating k = ", k, " in 0:", hSystem.size)
    p = mod(k + q, hSystem.size)

    ### do it?
        input = OrderedDict(
            "system size" => n,
            "momentum sector" => p,
            "magnetization sector" => 0,
            "coupling constant" => J,
            "magnon interaction" => β,
            "intrachain coupling" => K,
            "intrachain isotropy" => α,
            "hopping constant" => t
        )

        @time tJSystem, tJBasis, tJModel, tJFactorization = Main.SingleHole.run(input, howmany = 2, factor = true)

        vals, vecs, _ = tJFactorization
        mc = getMc(vecs[1], tJBasis, tJSystem)

        ### !now we may have 2-fold degeneracy (remove spin up or down)
        epsilon = 10^-12
        ΔE = vals[2] - vals[1]
        tMc = zeros(Float64, 2*tJSystem.size - 1)
        if abs(ΔE) < epsilon
            tMc = getMc(vecs[2], tJBasis, tJSystem)
            mc = 0.5 .* (mc .+ tMc)
        end

        println()

    return MC(n, β, J, t, k, K, α, collect(-(tJSystem.size - 1) : (tJSystem.size - 1)), mc)
end

### TODO: sum over magnetic configurations of other chains
function getMc(stateVector::Vector{ComplexF64}, basis, system)
    result = zeros(Float64, 2*system.size - 1)
    for m in -(system.size - 1) : (system.size - 1)
        state = zeros(Int64, system.size^2)
        if m < 0
            for it in (system.size + m):(system.size - 1)
                state[system.size^2 - it] = 1
            end
        else
            for it in 1:m
                state[system.size^2 - system.size + it] = 1
            end
        end
        if mod(m, 2) == 1
            state[system.size^2] = 1
        end
        binaryState = sum([state[it] * 2^(it-1) for it in 1:length(state)])
        result[system.size + m] = abs2(stateVector[basis[binaryState]])
    end
    ### rotate sublattice
    ### sum over 0 magnetization subspace
    ### rotate back

    return result
end

function getMcAll(stateVector::Vector{ComplexF64}, basis, system)
    result = zeros(Float64, 2*system.size - 1)
    for (state, _) in basis
        mask = Main.SingleHole.getMask(system)
        spinState = Main.SingleHole.sublatticeRotation(state, mask)
        if spinState < 2^(system.size^2 - system.size + 1)
            for m in -(system.size - 1) : (system.size - 1)
                # state = zeros(Int64, system.size^2)
                if m < 0
                    for it in (system.size + m):(system.size - 1)
                        state[system.size^2 - it] = 1
                    end
                else
                    for it in 1:m
                        state[system.size^2 - system.size + it] = 1
                    end
                end
                if mod(m, 2) == 1
                    state[system.size^2] = 1
                end
                binaryState = sum([state[it] * 2^(it-1) for it in 1:length(state)])
                result[system.size + m] += abs2(stateVector[basis[binaryState]])
            end
        end
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
            "intrachain" => data[it].intrachain,
            "isotropy" => data[it].isotropy,
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
n = 4
β = 1.0
KRange = [1.0, 0.1, 0.01, 10^-10] .* J
αRange = [0.0]

if length(ARGS) > 0
    t = eval(Meta.parse(ARGS[1]))
end
if length(ARGS) > 1
    J = eval(Meta.parse(ARGS[2]))
end
if length(ARGS) > 2
    KRange = [n for n in eval(Meta.parse(ARGS[3]))]
end
if length(ARGS) > 3
    αRange = [x for x in eval(Meta.parse(ARGS[4]))]
end

parameters = Vector{Parameters}()
for K in KRange
    for α in αRange
        push!(parameters, Parameters(n, β, J, t, K, α))
    end
end

data = Vector{MC}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_mc(parameters[it])
end

saveData(data)
