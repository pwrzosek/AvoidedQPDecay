using OrderedCollections
using JSON

# if !any(names(Main, imported = true) .== :Heisenberg)
    include("heisenberg.jl")
# end

n = 16
J = 1.0
β = 1.0
μ = 1 + div(n, 2)

E = Inf
V = []
k = 0
for momentum in 0:(n-1)
    input = OrderedDict(
        "system size" => n,
        "momentum sector" => momentum,
        "magnetization sector" => μ,
        "coupling constant" => J,
        "magnon interaction" => β
    )
    file = open("input.json", "w")
    JSON.print(file, input, 1)
    close(file)
    print(">>> ")
    @time system, fc = Main.Heisenberg.run()
    vals, vecs, info = fc
    if info.converged == 0
        println("Warning! Calculation did not converge.")
        println()
    end
    println("Energy: ", vals[1])
    if vals[1] < E
        global E = vals[1]
        global v = vecs[1]
        global k = 2 * momentum / n
    end
end
println("\n", "Ground state found for momentum: ", k, " π")
println("Ground state energy: ", E)


# @time sys, fc = Main.Heisenberg.run()

### !Important note:
### for system.size = 4l        -> GS @ k = 0
### for system.size = 4l + 2    -> GS @ k = π
