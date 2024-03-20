outputDirectory = "./reports/"
inputFileName = "qpd.sh"
nThreads = 2

s = 3
L = 400
for j in [0.3 0.25 0.2 0.15 0.085 0.075 0.065]
    for tp in [-0.3 0.0 0.3]
        outputFileName = string("qpd-s=", s, "-L=", L, "-j=", j, "-tp=", tp, ".out")
	outputCmd = string("#SBATCH --output=\"", outputDirectory, outputFileName, "\"")

        file = open(inputFileName, "w")

        println(file, "#!/bin/bash -l")
        println(file, "#SBATCH -J qpd-jl")
        println(file, "#SBATCH -N 1")
        println(file, "#SBATCH --ntasks-per-node 1")
        println(file, "#SBATCH --mem 8000")
        println(file, "#SBATCH -A G73-29")
        println(file, "#SBATCH -p topola")
        println(file, outputCmd)
        println(file)
        println(file, "module load apps/julia/1.1.1")
        println(file, "export JULIA_NUM_THREADS=", nThreads)
        println(file)
        println(file, "srun -n 1 julia main.jl s= ", s, " L= ", L, " j= ", j, " tp= ", tp)

        close(file)

        run(`sbatch qpd.sh`, wait = true)
    end
end
