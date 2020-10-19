outputDirectory = "../../reports/"
inputFileName = "qpd.sh"
nThreads = 12

outputFileName = string("1DtJ.out")
outputCmd = string("#SBATCH --output=\"", outputDirectory, outputFileName, "\"")

file = open(inputFileName, "w")

println(file, "#!/bin/bash -l")
println(file, "#SBATCH -J qpd-jl")
println(file, "#SBATCH -N 1")
println(file, "#SBATCH --ntasks-per-node 1")
println(file, "#SBATCH --mem 60000")
println(file, "#SBATCH -A G73-29")
println(file, "#SBATCH -p topola")
println(file, outputCmd)
println(file)
println(file, "module load apps/julia/1.1.1")
println(file, "export JULIA_NUM_THREADS=", nThreads)
println(file)
println(file, "srun -n 1 julia main.jl")

close(file)

run(`sbatch qpd.sh`, wait = true)
