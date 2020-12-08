# if !any(names(Main, imported = true) .== :Heisenberg)
    include("heisenberg.jl")
# end

sys, bs = @time Main.Heisenberg.run()
