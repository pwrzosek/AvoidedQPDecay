include("setup.jl")
include("half_filling.jl")

input = OrderedDict(
    "system size" => 20,
    "momentum sector" => 0,
    "magnetization sector" => 0,
    "coupling constant" => 1.0,
    "magnon interaction" => 1.0,
    "hopping constant" => 1.0
)

@time hf_run(input, howmany = 1)
