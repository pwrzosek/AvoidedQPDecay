include("declarations.jl")

# deltaRange = [0.0001]
# tPrimeRange = [it for it in -0.3:0.3:0.3]
# jRange = [0.4, 0.1, 0.05]
#
# for itd in 1:length(deltaRange)
#     global g_iDelta = deltaRange[itd] * im
#     z = zeros(length(jRange), length(tPrimeRange))
#     for itj in 1:length(jRange)
#         global g_j = jRange[itj]
#         for itt in 1:length(tPrimeRange)
#             global g_tPrime = g_t * tPrimeRange[itt]

            initialize()

            # k-dependent
            getGreensFunctionTemplate(false)

            # local
            # calculateLocalSpectralFunctionTemplate()

            finitialize()

            # drawLocalSpectralFunction(0.001, -4, 8, 2)
            # saveFigure(3, 2, "0")
            # writeLocalSpectralFunction(string("local", g_numberOfSites, "_", Int64(10 * g_j), ".txt"), 0.001, -3, 8)

            # if g_j < 0.05
            #     directory = string("J=", g_j, "_t'=", g_tPrime)
            #     draw() # (case = 0, kPath = g_kPath)
            #     saveFigure(0, 2, "0", directory)
            #     for (i, k) in enumerate(g_kPath)
            #            draw(0, k)
            #            saveFigure(i, 2, "0", directory)
            #     end
            # end

    #         z[itj, itt] = sqrt(getGroundStateResidue()[1][1])
    #         # draw()
    #     end
    # end

    # labels = reshape([string("t' = ", tPrimeRange[jt]) for jt in 1:length(tPrimeRange)], (1, length(tPrimeRange)))
    # display(Plots.plot(jRange, z, xlim = (-0.001, 0.031), ylim = (-0.001, 0.021), markershape = :auto, framestyle = :box, label = labels, xlabel = "J / t", ylabel = "z(t', J)", title = string("t-t'-Jz"), dpi = 600))
    # outputFileName = string("z0", "_L=", g_maxKrylovSpaceDimension)
    # saveFigure(0, 2, "z", outputFileName, false, false)

#     directory = string("z0")
#     if !isdir(directory)
#         mkdir(directory)
#     end
#     fileName = string(directory, "/", g_numberOfSites, "_L=", g_maxKrylovSpaceDimension, "_d=", imag(g_iDelta), ".txt")
#     file = open(fileName, "w")
#     writedlm(file, z, '\t')
#     close(file)
# end
