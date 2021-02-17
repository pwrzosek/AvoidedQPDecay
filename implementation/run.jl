using OrderedCollections
using LinearAlgebra
using JSON
using Dates
using Printf
using Plots

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
q = 0 + div(n, 2) * div(mod(n, 4), 2)

ωRange = collect(J/2-3:0.005:J/2+7)
kRange = (2 / n) .* collect(0:n)
βRange = collect(1.0:-0.1:1.0)
frames = Vector{Any}(undef, length(βRange))

for it in 1:length(βRange)
β = βRange[it]

### Heisenberg Ground State
input = OrderedDict(
    "system size" => n,
    "momentum sector" => q,
    "magnetization sector" => 1 + div(n, 2),
    "coupling constant" => J,
    "magnon interaction" => β
)
# file = open(dir["heisenberg"] * "input.json", "w")
# JSON.print(file, input, 1)
# close(file)

println()
@time hSystem, hBasis, hFactorization = Main.Heisenberg.run(input)

vals, vecs, info = hFactorization
GSE, GSV = vals[1], vecs[1]

println("\n", info)

### Spectral Function
spectrum = Array{Float64, 2}(undef, length(ωRange), hSystem.size + 1)
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
    # file = open(dir["tJModel"] * "input.json", "w")
    # JSON.print(file, input, 1)
    # close(file)

    @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input)

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
                coefficient = GSV[coordinate] * phase * sqrt(periodicity / l^2)
                initialState[tJBasis[repState]] += coefficient
            end
        end
    end
    initialState

    ### lanczos method for spectral function
    @time spectrum[:, k + 1] .= Main.SpectralFunction.run(ωRange .+ GSE, 0.05im, initialState, tJModel)
    println()
end

println("-- DONE @ ", string(Time(now()), " (", today(), ")"), " --\n")

frames[it] = spectrum

end

function heatgif(A::AbstractArray, kRange, ωRange, βRange; kwargs...)
    p = heatmap(kRange, ωRange, A[1])
    anim = @animate for it in append!(collect(1:length(A)), collect(length(A):-1:1))
        println("> Plotting ", it, " out of ", length(A), ".")
        heatmap!(p[1], kRange, ωRange, A[it], title = string("β = ", @sprintf "%.2f" βRange[it]); kwargs...)
    end
    return anim
end

# anim = heatgif(frames, kRange, ωRange, βRange, clim = (0, 1), colorbar = true, xlabel = "k / π", ylabel = "ω / t")
# gif(anim, "mmint_20.gif", fps = 8)

h = [heatmap(kRange, ωRange, 2*frames[it], clim = (0, 1), colorbar = true, title = string("β = ", βRange[it]), xlabel = "k / π", ylabel = "ω / t")
    for it in 1:length(frames)]

# for it in 1:length(h)
#     name = "plots/spectra/" * string(n) * "_s_β=" * string(βRange[it]) * ".png"
#     png(h[it], name)
# end
for it in 1:length(h)
    display(h[it])
end

# [sum(frames[1][:, k+1])*0.005 for k in 0:n]
