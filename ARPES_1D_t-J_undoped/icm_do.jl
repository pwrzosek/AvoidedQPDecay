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
include(dir["tJModel"] * "dos.jl")

struct DOS
    size::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
    momentum::Vector{Float64}
    spectrum::Array{Float64, 2}
    weights::Array{Float64, 2}
end

struct Parameters
    n::Int64
    β::Float64
    J::Float64
    t::Float64
end

function run_d(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    β = parameters.β
    J = parameters.J
    t = parameters.t

    howmany = 20

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


    ### DOS
    tJEigenvalues = Array{Float64, 2}(undef, howmany, hSystem.size + 1)
    weights = Array{Float64, 2}(undef, howmany, hSystem.size + 1)
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

        @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input, howmany = howmany)

        initialState = getInitialState(GSV, hBasis, hSystem, tJBasis, tJSystem)

        vals, vecs, _ = tJFactorization
        tJEigenvalues[:, k + 1] = vals[1:howmany] .- GSE

        ### calculate overlaps with HGS
        for it in 1:howmany
            weights[it, k + 1] = abs(dot(initialState, vecs[it]))^2
        end

        println()
    end

    return DOS(n, β, J, t, kRange, tJEigenvalues, weights) #DOS(n, β, J, t, kRange, ωRange, dos)
end

function saveData(data::Vector{DOS})
    dosData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Float64}, Array{Float64, 2}}}}(undef, length(data))

    for it in 1:length(data)
        dosData[it] = OrderedDict(
            "size" => data[it].size,
            "interaction" => data[it].interaction,
            "coupling" => data[it].coupling,
            "hopping" => data[it].hopping,
            "momentum" => data[it].momentum,
            "spectrum" => data[it].spectrum,
            "weights" => data[it].weights
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/dos", tail, ".json"), "w")
    JSON.print(file, dosData, 1)
    close(file)

    return nothing
end

function getInitialState(GSV, hBasis, hSystem, tJBasis, tJSystem)::Vector{Complex{Float64}}
    l::Int = hSystem.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    q = hSystem.momentum
    p = tJSystem.momentum

    ### not proved
    k = mod(p - q, hSystem.size)

    ### spin of removed electron: 1 -> up, 0 -> down
    removedSpin = 0

    ### set sublattice rotation masks
    mask = sum(1 << site for site in 0:2:(hSystem.size - 1))

    ### calculate initial state
    initialState = zeros(Complex{Float64}, length(tJBasis))
    for (state, coordinate) in hBasis
        spinState = xor(state, mask)
        for position in 1:l
            R = position - 1
            RValue = 1 << R
            repState = state
            for _ in 1:position
                repState = Main.Heisenberg.bitmov(repState, l, false, hb = highestBit, hv = highestValue)
            end
            phase = exp(2π * im * q * R / l) * exp(2π * im * p / l)
            periodicity = Main.Heisenberg.getPeriodicity(repState, hSystem)
            coefficient = sqrt(0.5) * GSV[coordinate] * phase * sqrt(periodicity / l^2)
            initialState[tJBasis[repState]] += coefficient
        end
    end
    return initialState
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

data = Vector{DOS}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_d(parameters[it])
end

saveData(data)
