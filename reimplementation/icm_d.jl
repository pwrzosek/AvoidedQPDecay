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
    energy::Vector{Float64}
    spectrum::Array{Float64, 2}
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

    vals, _, info = hFactorization
    GSE = vals[1]

    println("\n", info)

    ### Spectral Function
    dos = Array{Float64, 2}(undef, length(ωRange), hSystem.size + 1)
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

        @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input)

        tJEigenvalues, _, _ = tJFactorization

        ### lanczos method for spectral function
        @time dos[:, k + 1] .= Main.Dos.run(ωRange .+ GSE, iDelta, tJEigenvalues)

        println()
    end

    return DOS(n, β, J, t, kRange, ωRange, dos)
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
            "energy" => data[it].energy,
            "spectrum" => data[it].spectrum
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("./data/dos", tail, ".json"), "w")
    JSON.print(file, dosData, 1)
    close(file)

    return nothing
end

### System Parameters
t = 1.0
J = 0.4
nRange = [n for n in 8:2:8]
βRange = [-1.0, 0.0, 1.0]

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
