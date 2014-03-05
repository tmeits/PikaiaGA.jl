
module Pikaia

importall Base

export 
    pikaia, 
    rqsort

function rqsort(n:: Int, a:: Array{Float64, 1}, p:: Array{Int, 1})
# Return integer array p which indexes array a in increasing order.
# Array a is not disturbed.
    p=sortperm(a) # the permutation to sort an array
return p
end

function urand()
# Return the next pseudo-random deviate from a sequence which is
# uniformly distributed in the interval [0,1]
return rand()
end

function setctl(ctrl:: Array{Float64, 1}, n:: Int)
# Set control variables and flags from input and defaults
const DEFAULT = [100, 500, 5, .85, 2, .005, .0005, .25, 1, 1, 1, 0]
    for i=1:length(ctrl)
        if ctrl[i] < 0.
            ctrl[i] = DEFAULT[i]
        end
    end
    np     = int(ctrl[1])
    ngen   = int(ctrl[2])
    nd     = int(ctrl[3])
    pcross = ctrl[4]
    imut   = int(ctrl[5])
    pmut   = ctrl[6]
    pmutmn = ctrl[7]
    pmutmx = ctrl[8]
    fdif   = ctrl[9]
    irep   = int(ctrl[10])
    ielite = int(ctrl[11])
    ivrb   = int(ctrl[12])
    status = 0

# Print a header
    if ivrb > 0
        @printf("******************************************************************\n")
        @printf("*            Pikaia Genetic Algorithm Report                     *\n")
        @printf("******************************************************************\n")
        @printf("   Number of Generations evolving: %4i\n", ngen)
        @printf("       Individuals per generation: %4i\n", np)
        @printf("    Number of Chromosome segments: %4i\n", n)
        @printf("    Length of Chromosome segments: %4i\n", nd)
        @printf("            Crossover probability: %9.4f\n", pcross)
        @printf("            Initial mutation rate: %9.4f\n", pmut)
        @printf("            Minimum mutation rate: %9.4f\n", pmutmn)
        @printf("            Maximum mutation rate: %9.4f\n", pmutmx)
        @printf("    Relative fitness differential: %9.4f\n", fdif)
       
        if imut == 1
            @printf("                    Mutation Mode: Uniform, Constant Rate\n")
        elseif imut == 2
            @printf("                    Mutation Mode: Uniform, Variable Rate (F)\n")
        elseif imut == 3
            @printf("                    Mutation Mode: Uniform, Variable Rate (D)\n")
        elseif imut == 4
            @printf("                    Mutation Mode: Uniform+Creep, Constant Rate\n")
        elseif imut == 5
            @printf("                    Mutation Mode: Uniform+Creep, Variable Rate (F)\n")
        elseif imut == 6
            @printf("                    Mutation Mode: Uniform+Creep, Variable Rate (D)\n")
        end
        if irep == 1 
            @printf("                Reproduction Plan: Full generational replacement\n")
        elseif irep == 2
            @printf("                Reproduction Plan: Steady-state-replace-random\n")
        elseif irep == 3
            @printf("                Reproduction Plan: Steady-state-replace-worst\n")
        end
    end
# Check some control values
    if imut != 1 && imut != 2 && imut != 3 && imut != 4 && imut != 5 && imut != 6
        @printf(" ERROR: illegal value for imut (ctrl(5))\n")
        status = 5
    end

    if fdif > 1
        @printf(" ERROR: illegal value for fdif (ctrl(9))\n")
        status = 9
    end

    if irep != 1 && irep != 2 && irep != 3
        @printf(" ERROR: illegal value for irep (ctrl(10))\n")
        status = 10
    end

    if pcross > 1.0 || pcross < 0.0
        @printf(" ERROR: illegal value for pcross (ctrl(4))\n")
        status = 4
    end

    if ielite  != 1.0 && ielite != 1
        @printf(" ERROR: illegal value for ielite (ctrl(11))\n")
        status = 11
    end

    if irep == 1 && imut == 1 && pmut > 0.5 && ielite == 0
        @printf(" WARNING: dangerously high value for pmut (ctrl(6));\n")
        @printf(" (Should enforce elitism with ctrl(11)=1.)\n")
    end

    if irep == 1 && imut == 2 && pmutmx > 0.5 && ielite == 0
        @printf(" WARNING: dangerously high value for pmutmx (ctrl(8));\n")
        @printf(" (Should enforce elitism with ctrl(11)=1.)\n")
    end

    if fdif < 0.33 && irep != 3
        @printf(" WARNING: dangerously low value of fdif (ctrl(9))\n") 
    end

    if mod(np, 2) > 0
        np = np - 1
        @printf("WARNING: decreasing population size (ctrl(1)) to np= %4i\n", np)
    end

    return (status, np, ngen, nd, imut, irep, ielite, ivrb, 
        pcross, pmutmn, pmutmx, pmut, fdif)

end

function  rnkpop(n:: Int, arrin:: Array{Float64, 1})
# Compute the key index
    rank = rand(n)
    indx = rqsort(n, arrin, [1: n])
# ...and the rank order
    for i = 1 : n
        rank[indx[i]] = n - i + 1
    end
    return (indx, rank)
end

function pikaia(ff:: Function, n:: Int, ctrl:: Array{Float64, 1})
# Optimization (maximization) of user-supplied "fitness" function
# ff  over n-dimensional parameter space  x  using a basic genetic
# algorithm method.

# Version 0.0.1   [ 2014 February 21 ]

const NMAX = 32   # NMAX is the maximum number of adjustable parameters (n <= NMAX)
const PMAX = 128  # PMAX is the maximum population (ctrl(1) <= PMAX)
const DMAX = 6    # DMAX is the maximum number of Genes (digits) per Chromosome segement (parameter) (ctrl(3) <= DMAX)

#  Local variables


# Init output:
    x = Array{Float64, n}
    x = rand(n)

    f :: Float64 = 0.
    status :: Int = 0

# Set control variables from input and defaults
    (status, np, ngen, nd, imut, irep, ielite, ivrb, 
        pcross, pmutmn, pmutmx, pmut, fdif) = setctl(ctrl, n)
    if status != 0
        println(" Control vector (ctrl) argument(s) invalid")
        return (x, f, status)
    end
# Make sure locally-dimensioned arrays are big enough
    if  n > NMAX || np > PMAX || nd > DMAX 
        println(" Number of parameters, population, or genes too large")
        status = -1
        return (x, f, status)
    end

# Compute initial (random but bounded) phenotypes
    for ip = 1 : np
        for k = 1 : n
            oldph[k, ip] = urand()
        end
    end
# Rank initial population by fitness order
    (ifit, jfit) = rnkpop(np, fitns)

# Main Generation Loop
    for ig = 1 : ngen
# Main Population Loop
        newtot = 0
        for ip = 1 : np / 2
        end # End of Main Population Loop
    end # End of Main Generation Loop 

# Return in typle
return (x, f, status)
end # pikaia

end # Pikaia
#
# #3F3434 #DFE8E2
# http://habrahabr.ru/post/161009/
# https://github.com/johnmyleswhite/HopfieldNets.jl
# http://habrahabr.ru/post/125999/
