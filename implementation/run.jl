using OrderedCollections
using LinearAlgebra
using JSON
using Dates

println(">> Today: ", today(), "   Time: ", Time(now()), " <<\n")

dir = Dict(
    "heisenberg" => "heisenberg/code/julia/",
    "tJModel" => "tJ_single_hole/code/julia/"
)

include(dir["heisenberg"] * "heisenberg.jl")
include(dir["tJModel"] * "tJmodel.jl")
include(dir["tJModel"] * "spectral.jl")

### System Parameters
n = 16
t = 1.0
J = 0.4
β = 1.0
q = 0

### Heisenberg Ground State
input = OrderedDict(
    "system size" => n,
    "momentum sector" => q,
    "magnetization sector" => 1 + div(n, 2),
    "coupling constant" => J,
    "magnon interaction" => β
)
file = open(dir["heisenberg"] * "input.json", "w")
JSON.print(file, input, 1)
close(file)

println()
@time hSystem, hBasis, hFactorization = Main.Heisenberg.run()

vals, vecs, info = hFactorization
GSE, GSV = vals[1], vecs[1]

println("\n", info)

### Spectral Function
test = []
for k in 0:hSystem.size
    println("Evaluating k = ", k, " in 0:", hSystem.size)
    p = mod(k + q, hSystem.size)

    input = OrderedDict(
        "system size" => n,
        "momentum sector" => p,
        "magnetization sector" => 1 + div(n, 2),
        "coupling constant" => J,
        "magnon interaction" => β,
        "hopping constant" => t
    )
    file = open(dir["tJModel"] * "input.json", "w")
    JSON.print(file, input, 1)
    close(file)

    @time tJSystem, tJBasis, tJModel = Main.tJmodel.run()

    ### calculate initial state
    l::Int = hSystem.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1
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
                coefficient = GSV[coordinate] * phase * sqrt(periodicity / l)
                initialState[tJBasis[repState]] += coefficient
            end
        end
    end

    ### lanczos method for spectral function
    ωRange = collect(-4:0.01:6)
    @time Main.SpectralFunction.run(ωRange, initialState, tJModel)
    println()
end

println("-- DONE @ ", string(Time(now()), " (", today(), ")"), " --\n")
