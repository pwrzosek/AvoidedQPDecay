using OrderedCollections
using KrylovKit
using LinearMaps

struct System
    size::Int64
    momentum::Int64
    magnetization::Int64
    coupling::Float64
    interaction::Float64
    hopping::Float64
end

function make_system(input::OrderedDict)::System
    system::System = System(
        input["system size"],
        input["momentum sector"],
        input["magnetization sector"],
        input["coupling constant"],
        input["magnon interaction"],
        input["hopping constant"]
    )
    return system
end

Basis = OrderedDict{Int64,Int64}
