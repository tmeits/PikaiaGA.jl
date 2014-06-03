# include("testfunctions.jl")
# include("pikaia.jl")

module TestFunctions
# Test functions for optimization
# http://en.wikipedia.org/wiki/Test_functions_for_optimization

export 
    rastrigin,
    rastriginsfcn

# *********************************************************************
function rastriginsfcn(x::Vector{Float64})
# =====================================================================
# Rastrigin function of two variables
# =====================================================================

    if length(x) == 2
        20+x[1]^2+x[2]^2-10*(cos(2*pi*x[1])+cos(2*pi*x[2]))
    else
        error("number of input arguments id bad")
    end

end

function rastrigin(x::Vector{Float64})
# http://en.wikipedia.org/wiki/Rastrigin_function

    a = 10.0
    n = length(x)

    fitness = a*n
    for i = 1:n
        fitness += x[i]^2 - a*cos(2*pi*x[i])
    end

    return fitness
end

end
# https://code.google.com/p/psomatlab/
# http://jenyay.net/Programming/ParticleSwarm
