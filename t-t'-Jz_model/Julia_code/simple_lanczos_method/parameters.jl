g_latticeType = 4

g_t = 1.
g_j = g_t * 0.4
g_tPrime = g_t * 0.

g_iDelta = 0.05 * im

g_maxKrylovSpaceDimension = 200

if length(ARGS) > 0
    for (i, arg) in enumerate(ARGS)
        if arg == "s="
            global g_latticeType = parse(Int64, ARGS[i+1])
        end
        if arg == "tp="
            global g_tPrime = parse(Float64, ARGS[i+1])
        end
        if arg == "j="
            global g_j = parse(Float64, ARGS[i+1])
        end
        if arg == "L="
            global g_maxKrylovSpaceDimension = parse(Int64, ARGS[i+1])
        end
        if arg == "d="
            global g_iDelta = parse(Float64, ARGS[i+1]) * im
        end
    end
end

nothing

# example:
#   user$ julia main.jl s= 2 j= 0.05
