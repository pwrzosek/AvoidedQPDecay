using Dates
using Plots
using Plots.PlotMeasures
using Printf
using LinearAlgebra
using StaticArrays
using DelimitedFiles

BLAS.set_num_threads(1)

g_showMessages = false
g_showProgressBar = false

g_isParallel = false

if g_isParallel # warning: gc causes problems when run on ICM machine
    g_showMessages = false
    g_showProgressBar =  false
end

nothing
