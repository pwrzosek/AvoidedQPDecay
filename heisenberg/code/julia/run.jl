# if !any(names(Main, imported = true) .== :Heisenberg)
    include("heisenberg.jl")
# end

basis = @time Main.Heisenberg.run()
