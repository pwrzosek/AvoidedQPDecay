using Dates
using Plots
using Plots.PlotMeasures
using Printf
using LinearAlgebra
using StaticArrays
using DelimitedFiles

g_showMessages = true
g_showProgressBar = true

g_isParallel = false

if g_isParallel # warning: gc causes problems when run on ICM machine
    g_showMessages = false
    g_showProgressBar =  false
end

nothing
