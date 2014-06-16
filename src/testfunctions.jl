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

# Одномерная функция тестирования
function f(x)

    return (abs(x) + cos(x))
end

# Min x = -20
# Max x = 20
# x = 0 -> 1

function Haupt_and_Haupt_2004(x)

#rescale input variables:
     

    return (x^2 + x) * cos(x)
end    

# Min x = -10
# Max x = 10
# Fitness function value = 47.70562
# Solution x =  6.560539

# http://www.jstatsoft.org/v53/i04/paper


# *********************************************************************
function twod(x)
# =====================================================================
# Compute sample fitness function (2-d landscape)
# =====================================================================

#     Constant
#    const pi     =
    const sigma2 = 0.15
    const nn     = 9

    rr = 0.0         

    if x[1] > 1.0 || x[2] > 1.0 
        error("stop")
    end

    rr = sqrt( (0.5-x[1]) ^ 2 + (0.5-x[2]) ^ 2 )
    
    return cos(rr*nn*pi) ^ 2 * exp((-rr ^ 2) / sigma2)
end


