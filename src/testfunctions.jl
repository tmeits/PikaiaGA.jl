#
module TestFunctions

# *********************************************************************
function rastriginsfcn(x::Vector{Float64})
# =====================================================================
# Rastrigin's function, a function that is often used to test the genetic algorithm.
# http://www.mathworks.com/help/gads/rastrigins_with_line2.png
# http://www.mathworks.com/help/gads/ras_contour.png
# http://www.mathworks.com/help/gads/eqn1198784102.png
# =====================================================================

    if length(x) == 2
        20+x[1]^2+x[2]^2-10*(cos(2*pi*x[1])+cos(2*pi*x[2]))
    else
        error("number of input arguments id bad")
    end

end

end
# https://code.google.com/p/psomatlab/
# http://jenyay.net/Programming/ParticleSwarm
