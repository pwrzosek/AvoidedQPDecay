include("single_hole.jl")

@time system, basis, model, factorization = Main.SingleHole.run()
