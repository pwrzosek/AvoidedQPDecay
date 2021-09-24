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
    "tJModel" => "tJ_single_hole/code/julia/",
)

include(dir["heisenberg"] * "heisenberg.jl")
include(dir["tJModel"] * "tJmodel.jl")
include(dir["tJModel"] * "spectral.jl")

struct QP
    energy::Float64
    weight::Float64
    gap::Float64
    size::Int64
    momentum::Int64
    interaction::Float64
    coupling::Float64
    hopping::Float64
end

struct Parameters
    n::Int64
    k::Int64
    β::Float64
    J::Float64
    t::Float64
end

function run_z(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    k = parameters.k
    β = parameters.β
    J = parameters.J
    t = parameters.t
    q = 0 + div(n, 2) * div(mod(n, 4), 2)
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => q,
        "magnetization sector" => 1 + div(n, 2),
        "coupling constant" => J,
        "magnon interaction" => β
    )

    ### Calculate Heisenberg GS
    println()
    @time hSystem, hBasis, hFactorization = Main.Heisenberg.run(input)

    vals, vecs, info = hFactorization

    GSE, GSV = vals[1], vecs[1]

    ### Spectrum and z

    ### t-J model input parameters
    p = mod(k + q, hSystem.size)
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => p,
        "magnetization sector" => 1 + div(n, 2),
        "coupling constant" => J,
        "magnon interaction" => β,
        "hopping constant" => t
    )

    tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input)

    initialState = getInitialState(GSV, hBasis, hSystem, tJBasis, tJSystem)

    vals, vecs, info = tJFactorization

    tJGSE, tJGSV = vals[1], vecs[1]

    pole = tJGSE - GSE
    residue = abs(dot(tJGSV, initialState))^2

    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)
    @time spectrum = Main.SpectralFunction.run(ωRange .+ GSE, iDelta, initialState, tJModel)

    return (QP(pole[1], residue[1], vals[2] - vals[1], n, k, β, J, t), ωRange, spectrum, tJGSV)
end

function getInitialState(GSV, hBasis, hSystem, tJBasis, tJSystem)
    l::Int = hSystem.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    q = hSystem.momentum
    p = tJSystem.momentum

    ### calculate initial state
    initialState = zeros(Complex{Float64}, length(tJBasis))
    for (state, coordinate) in hBasis
        for position in 1:l
            R = position - 1
            RValue = 1 << R
            if (state & RValue) == 0 ## == 0 for spin down annihilation or != for spin up (second one not implemented yet)
                repState = state
                for _ in 1:position
                    repState = Main.Heisenberg.bitmov(repState, l, false, hb = highestBit, hv = highestValue)
                end
                phase = exp(2π * im * q * R / l) * exp(2π * im * p / l)
                periodicity = Main.Heisenberg.getPeriodicity(repState, hSystem)
                coefficient = GSV[coordinate] * phase * sqrt(periodicity / l^2)
                initialState[tJBasis[repState]] += coefficient
            end
        end
    end
    return initialState
end

function saveData(data::Vector{Tuple{QP, Vector{Float64}, Vector{Float64}, Vector{Complex{Float64}}}})
    qpData = Vector{OrderedDict{String, Union{Int64, Float64}}}(undef, length(data))
    sData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Vector{Float64}}}}}(undef, length(data))
    gsData = Vector{OrderedDict{String, Union{Int64, Float64}}}(undef, length(data))

    # calculate overlaps
    overlaps = Vector{Float64}(undef, length(data))
    referenceIndices = findall(x -> x[1].interaction == 1.0, data)

    doOverlaps = false
    if length(referenceIndices) > 0
        doOverlaps = true
    end

    if doOverlaps
        for it in referenceIndices
            matchingIndices = findall(x ->
                (x[1].size == data[it][1].size) &&
                (x[1].momentum == data[it][1].momentum) &&
                (x[1].coupling == data[it][1].coupling) &&
                (x[1].hopping == data[it][1].hopping),
            data)
            for jt in matchingIndices
                overlaps[jt] = abs(dot(data[it][4], data[jt][4]))^2
            end
        end
    end

    for it in 1:length(data)
        qpData[it] = OrderedDict(
            "energy" => data[it][1].energy,
            "weight" => data[it][1].weight,
            "gap" => data[it][1].gap,
            "size" => data[it][1].size,
            "momentum" => 2 * data[it][1].momentum / data[it][1].size,
            "interaction" => data[it][1].interaction,
            "coupling" => data[it][1].coupling,
            "hopping" => data[it][1].hopping
        )
        sData[it] = OrderedDict(
            "spectrum" => [data[it][2], data[it][3]],
            "size" => data[it][1].size,
            "momentum" => 2 * data[it][1].momentum / data[it][1].size,
            "interaction" => data[it][1].interaction,
            "coupling" => data[it][1].coupling,
            "hopping" => data[it][1].hopping
        )
        if doOverlaps
            gsData[it] = OrderedDict(
                "overlap" => overlaps[it],
                "size" => data[it][1].size,
                "momentum" => 2 * data[it][1].momentum / data[it][1].size,
                "interaction" => data[it][1].interaction,
                "coupling" => data[it][1].coupling,
                "hopping" => data[it][1].hopping
            )
        end
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/qp", tail, ".json"), "w")
    JSON.print(file, qpData, 1)
    close(file)

    file = open(string("./data/s", tail, ".json"), "w")
    JSON.print(file, sData, 1)
    close(file)

    if doOverlaps
        file = open(string("./data/gs", tail, ".json"), "w")
        JSON.print(file, gsData, 1)
        close(file)
    end

    return nothing
end

t = 1.0
J = 0.4
nRange = [n for n in 16:2:16]
βRange = [x for x in 0.0:1.0:1.0]

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

kDiv2pi = [0, 1/4, 1/2]
parameters = Vector{Parameters}()
for n in nRange
    kPoints = [0, div(n, 4), div(n, 2)]
    kAvail = isinteger.(n * kDiv2pi)
    for k in kPoints[kAvail]
        for β in βRange
            push!(parameters, Parameters(n, k, β, J, t))
        end
    end
end

data = Vector{Tuple{QP, Vector{Float64}, Vector{Float64}, Vector{Complex{Float64}}}}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_z(parameters[it])
end

saveData(data)
