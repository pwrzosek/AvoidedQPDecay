JRange = 0.0:0.01:1.0

for J in JRange
    try
        run(`julia main.jl s= 1 j= $J`, wait = true)
    catch
        J += 0.000000001
        run(`julia main.jl s= 1 j= $J`, wait = true)
    end
end
