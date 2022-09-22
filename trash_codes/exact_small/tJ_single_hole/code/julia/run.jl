include("tJmodel.jl")

# @time tJSystem, tJBasis, tJModel = Main.tJmodel.run()
@time sys, bs, md, fc = Main.tJmodel.run()
