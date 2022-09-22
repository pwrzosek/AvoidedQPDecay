using Test
using OrderedCollections
using Suppressor

include("single_hole.jl")

@suppress_err begin

minExactSolverSize = 4
maxExactSolverSize = 8
exactSolverStep = 4

minLanczosSolverSize = 8
maxLanczosSolverSize = 12
lanczosSolverStep = 4

println()

### Testing:
#   1.  GS is at momentum k = div(n, 4) --> YES
#   2.  GS is at magnetization m = 1 --> YES [m === abs(nSpinsUp - nSpinsDown)]
#   3.  GS is doubly degenerate --> YES
#   4.  GS is also at equivalent momentum point k = n - div(n, 4) --> YES
#   5.  GS at k = n - div(n, 4) is degenerate --> YES

@testset verbose = true "Heisenberg Exact Solver" begin
    @testset "size $i" for i in minExactSolverSize:exactSolverStep:maxExactSolverSize
        n = i
        GSE = Inf
        GSK = 0
        GSM = 0
        for k in 0:(n-1)
            for m in 1:2:(n-1)
                input = OrderedDict(
                    "system size" => n,
                    "magnetization sector" => m,
                    "momentum sector" => k,
                    "coupling constant" => 0.4,
                    "magnon fluctuations" => 1.0,
                    "magnon interaction" => 1.0,
                    "proximity interaction" => 1.0,
                    "hopping constant" => 1.0,
                    "solver type" => 1
                )
                _, _, _, factorization = Main.SingleHole.run(input, timed = false)
                if factorization[1][1] < GSE
                    GSE = factorization[1][1]
                    GSK = k
                    GSM = m
                end
            end
        end
        @test GSM == 1
        @test (GSK == div(n, 4)) || (GSK == n - div(n, 4))

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => n - div(n, 4),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 1.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 1
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE1 = factorization[1][1]

        @test factorization[1][1] ≈ factorization[1][2]

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => div(n, 4),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 1.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 1
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE2 = factorization[1][1]

        @test GSE1 ≈ GSE2
        @test factorization[1][1] ≈ factorization[1][2]
    end
end


### Testing:
#   1.  GS is at momentum k = div(n, 4) --> YES
#   2.  GS is at magnetization m = 1 --> YES [m === abs(nSpinsUp - nSpinsDown)]
#   3.  GS is doubly degenerate --> YES
#   4.  GS is also at equivalent momentum point k = n - div(n, 4) --> YES
#   5.  GS at k = n - div(n, 4) is degenerate --> YES

@testset verbose = true "Heisenberg Lanczos Solver" begin
    @testset "size $i" for i in minLanczosSolverSize:lanczosSolverStep:maxLanczosSolverSize
        n = i
        GSE = Inf
        GSK = 0
        GSM = 0
        for k in 0:(n-1)
            for m in 1:2:(n-1)
                input = OrderedDict(
                    "system size" => n,
                    "magnetization sector" => m,
                    "momentum sector" => k,
                    "coupling constant" => 0.4,
                    "magnon fluctuations" => 1.0,
                    "magnon interaction" => 1.0,
                    "proximity interaction" => 1.0,
                    "hopping constant" => 1.0,
                    "solver type" => 0
                )
                _, _, _, factorization = Main.SingleHole.run(input, timed = false)
                if factorization[1][1] < GSE
                    GSE = factorization[1][1]
                    GSK = k
                    GSM = m
                end
            end
        end
        @test GSM == 1
        @test (GSK == div(n, 4)) || (GSK == n - div(n, 4))

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => n - div(n, 4),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 1.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 0
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE1 = factorization[1][1]

        @test factorization[1][1] ≈ factorization[1][2]

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => div(n, 4),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 1.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 0
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE2 = factorization[1][1]

        @test GSE1 ≈ GSE2
        @test factorization[1][1] ≈ factorization[1][2]
    end
end


### Testing:
#   1.  GS is at momentum k = 0 --> YES
#   2.  GS is at magnetization m = 1 --> YES [m === abs(nSpinsUp - nSpinsDown)]
#   3.  GS is non-degenerate --> YES
#   4.  GS is also at momentum k = div(n, 2) --> YES
#   5.  GS at k = div(n, 2) is non-degenerate --> YES

@testset verbose = true "Ising Exact Solver" begin
    @testset "size $i" for i in minExactSolverSize:exactSolverStep:maxExactSolverSize
        n = i
        GSE = Inf
        GSK = 0
        GSM = 0
        for k in 0:(n-1)
            for m in 1:2:(n-1)
                input = OrderedDict(
                    "system size" => n,
                    "magnetization sector" => m,
                    "momentum sector" => k,
                    "coupling constant" => 0.4,
                    "magnon fluctuations" => 0.0,
                    "magnon interaction" => 1.0,
                    "proximity interaction" => 1.0,
                    "hopping constant" => 1.0,
                    "solver type" => 1
                )
                _, _, _, factorization = Main.SingleHole.run(input, timed = false)
                if factorization[1][1] < GSE
                    GSE = factorization[1][1]
                    GSK = k
                    GSM = m
                end
            end
        end
        @test GSM == 1
        @test GSK == 0 || GSK == div(n, 2)

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => div(n, 2),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 0.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 1
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE1 = factorization[1][1]

        @test abs(factorization[1][1] - factorization[1][2]) > 10^-12

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => 0,
            "coupling constant" => 0.4,
            "magnon fluctuations" => 0.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 1
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE2 = factorization[1][1]

        @test GSE1 ≈ GSE2
        @test abs(factorization[1][1] - factorization[1][2]) > 10^-12
    end
end


### Testing:
#   1.  GS is at momentum k = 0 --> YES
#   2.  GS is at magnetization m = 1 --> YES [m === abs(nSpinsUp - nSpinsDown)]
#   3.  GS is non-degenerate --> YES
#   4.  GS is also at momentum k = div(n, 2) --> YES
#   5.  GS at k = div(n, 2) is non-degenerate --> YES

@testset verbose = true "Ising Lanczos Solver" begin
    @testset "size $i" for i in minLanczosSolverSize:lanczosSolverStep:maxLanczosSolverSize
        n = i
        GSE = Inf
        GSK = 0
        GSM = 0
        for k in 0:(n-1)
            for m in 1:2:(n-1)
                input = OrderedDict(
                    "system size" => n,
                    "magnetization sector" => m,
                    "momentum sector" => k,
                    "coupling constant" => 0.4,
                    "magnon fluctuations" => 0.0,
                    "magnon interaction" => 1.0,
                    "proximity interaction" => 1.0,
                    "hopping constant" => 1.0,
                    "solver type" => 0
                )
                _, _, _, factorization = Main.SingleHole.run(input, timed = false)
                if factorization[1][1] < GSE
                    GSE = factorization[1][1]
                    GSK = k
                    GSM = m
                end
            end
        end
        @test GSM == 1
        @test GSK == 0 || GSK == div(n, 2)

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => div(n, 2),
            "coupling constant" => 0.4,
            "magnon fluctuations" => 0.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 0
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE1 = factorization[1][1]

        @test abs(factorization[1][1] - factorization[1][2]) > 10^-12

        input = OrderedDict(
            "system size" => n,
            "magnetization sector" => GSM,
            "momentum sector" => 0,
            "coupling constant" => 0.4,
            "magnon fluctuations" => 0.0,
            "magnon interaction" => 1.0,
            "proximity interaction" => 1.0,
            "hopping constant" => 1.0,
            "solver type" => 0
        )
        _, _, _, factorization = Main.SingleHole.run(input, timed = false)
        GSE2 = factorization[1][1]

        @test GSE1 ≈ GSE2
        @test abs(factorization[1][1] - factorization[1][2]) > 10^-12
    end
end

end
