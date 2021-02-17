using OrderedCollections
using LinearAlgebra
using Dates
using Plots

### upload needed libaries

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
    println("\n", ">> Today: ", today(), "   Time: ", Time(now()), " <<\n")

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
    println("> Heisenberg ground state calculation:")
    @time hSystem, hBasis, hFactorization = Main.Heisenberg.run(input)

    vals, vecs, info = hFactorization
    println("\n", info)

    GSE, GSV = vals[1], vecs[1]

    ### Spectrum and z
    println("> t-J model evaluation:")

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

    @time tJSystem, tJBasis, tJModel, tJFactorization = Main.tJmodel.run(input)

    initialState = getInitialState(GSV, hBasis, hSystem, tJBasis, tJSystem)

    vals, vecs, info = tJFactorization
    println("\n", info)

    tJGSE, tJGSV = vals[1], vecs[1]

    println("Norm initial : ", norm(initialState), " & norm GS : ", norm(tJGSV))

    pole = tJGSE - GSE
    residue = abs(dot(tJGSV, initialState))^2

    iDelta = 0.01im
    ωRange = collect(-3:0.001:7)
    spectrum = Main.SpectralFunction.run(ωRange .+ GSE, iDelta, initialState, tJModel)

    println("\n-- DONE @ ", string(Time(now()), " (", today(), ")"), " --\n")

    return (QP(pole[1], residue[1], n, k, β, J, t), ωRange, spectrum)
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

t = 1.0
J = 1.0
nRange = [4, 8, 12, 16, 20]
βRange = [1.0, 0.9, 0.5]

parameters = Vector{Parameters}()
for n in nRange
    kPoints = [div(n, 4)]
    for k in kPoints
        for β in βRange
            push!(parameters, Parameters(n, k, β, J, t))
        end
    end
end

data = Vector{Tuple{QP, Vector{Float64}, Vector{Float64}}}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_z(parameters[it])
end

function zPlots(data, kRange, nRange, βRange)
    nMax = nRange[end]
    plots = [plot() for _ in kRange]
    points = Array{Float64, 3}(undef, length(kRange), length(nRange), length(βRange))
    spectra = Array{Vector{Vector{Float64}}, 2}(undef, length(βRange), length(kRange))
    poles = Array{Float64, 2}(undef, length(βRange), length(kRange))
    for (qp, ωRange, spectrum) in data
        kIndex = findall(x -> x == 2π * qp.momentum / qp.size, kRange)[1]
        nIndex = findall(x -> x == qp.size, nRange)[1]
        βIndex = findall(x -> x == qp.interaction, βRange)[1]
        points[kIndex, nIndex, βIndex] = qp.weight
        if qp.size == nMax
            spectra[βIndex, kIndex] = [ωRange, spectrum]
            poles[βIndex, kIndex] = qp.energy
        end
    end
    for ik in 1:length(kRange)
        plots[ik] = plot(1 ./ nRange, points[ik, :, :],
            xlim = (0, 0.32), ylim = (0,1.0),
            labels = permutedims(string.(βRange)),
            frame = :box,
            xlabel = "1/L", ylabel = "z",
            markershape = :circle,
            markersize = 3.6,
            title = string("k = ", round(kRange[ik] / π, digits = 1), "π")
        )
        # display(plot(spectra[1, ik][1], spectra[1, ik][2], title = string(poles[1, ik])))
    end
    return plots, spectra
end

kRange = [π/2]
p, s = zPlots(data, kRange, nRange, βRange)

for ik in 1:length(kRange)
    display(p[ik])
    kName = string(round(kRange[ik] / π, digits = 1), "π")
    name = "plots/z/" * string(nRange[end]) * "_zgs_k=" * kName * ".png"
    png(p[ik], name)
end

# function ΔPlots(data, kRange, nRange, βRange)
#     nMax = nRange[end]
#     plots = [plot() for _ in kRange]
#     points = Array{Float64, 3}(undef, length(kRange), length(nRange), length(βRange))
#     for (qp, ωRange, spectrum) in data
#         kIndex = findall(x -> x == 2π * qp.momentum / qp.size, kRange)[1]
#         nIndex = findall(x -> x == qp.size, nRange)[1]
#         βIndex = findall(x -> x == qp.interaction, βRange)[1]
#         points[kIndex, nIndex, βIndex] = qp.energy
#     end
#     for ik in 1:length(kRange)
#         plots[ik] = plot(1 ./ nRange, points[ik, :, :], xlim = (0, 0.32), ylim = (0,1), labels = permutedims(string.(βRange)))
#         # display(plot(spectra[1, ik][1], spectra[1, ik][2], title = string(poles[1, ik])))
#     end
#     return plots
# end
#
# kRange = [0, π/2, π]
# p = ΔPlots(data, kRange, nRange, βRange)
#
# for ik in 1:length(kRange)
#     display(p[ik])
# end

# diagonal, offDiagonal, size, spectrum = Main.SpectralFunction.run(ωRange .+ GSE, iDelta, initialState, tJModel)
#
# f = []
# try
#     f = eigen(SymTridiagonal(diagonal, offDiagonal), 1:size)
# catch
#     f = eigen(SymTridiagonal(diagonal, offDiagonal), 1:size)
# end
# energy = f.values
# weight = abs(iDelta) * π * Main.SpectralFunction.spectralFunction(energy, iDelta, initialState, diagonal, offDiagonal, size)
#
# cutoff = 10^-10
# pole = energy[1]
# residue = weight[1]
# for it in 1:length(weight)
#     if weight[it] > cutoff
#         pole = energy[it]
#         residue = weight[it]
#         break
#     end
# end
