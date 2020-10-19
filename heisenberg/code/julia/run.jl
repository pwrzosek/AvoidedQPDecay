# if !any(names(Main, imported = true) .== :Heisenberg)
    include("heisenberg.jl")
# end

@time Main.Heisenberg.run()
