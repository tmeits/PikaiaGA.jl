
module Pikaia

importall Base

export 
    pikaia, 
    rqsort,
    init_pop,
    newpop,
    urand,
    encode,
    cross!,
    mutate!

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
st
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
# Calls external sort routine to produce key index and rank order
# of input array arrin (which is not altered).

#   Compute the key index
    rank = [1:n] #rand(n)
    indx = rqsort(n, arrin, [1: n])
#   ...and the rank order
    for i = 1 : n
        rank[indx[i]] = n - i + 1
    end
    return (indx, rank)
end

function select(np:: Int, jfit:: Array{Int, 1}, fdif:: Float64)
# Selects a parent from the population, using roulette wheel
# algorithm with the relative fitnesses of the phenotypes as
# the "hit" probabilities [see Davis 1991, chap. 1].
    
    idad = 0

    np1 = np + 1
    dice = urand()*np*np1
    rtfit = 0.
    for i= 1 : np
        rtfit = rtfit+np1+fdif*(np1-2*jfit[i])
        if rtfit >= dice
            idad=i
            break
        end
    end
    return idad
end

function init_pop(ff:: Function, n:: Int, np:: Int)
#   Compute initial (random but bounded) phenotypes
    old_ph = Array(Float64, n, np)
    fitns = Array(Float64, np)
#   old_ph = rand(n, np)
    for ip = 1 : np
        for k = 1 : n
            old_ph[k, ip] = urand()
        end
        fitns[ip] = ff(n, old_ph[:,ip])
    end
#   Rank initial population by fitness order
    (ifit, jfit) = rnkpop(np, fitns)
    return (old_ph, fitns, ifit, jfit)
end

function new_pop(ff:: Function, ielite:: Int, ndim:: Int, n:: Int, np:: Int, 
    oldph:: Array{Float64, 1}, fitns:: Array{Float64, 1})
# replaces old population by new; recomputes fitnesses & ranks
    
    nnew = np
    newph = oldph

#   if using elitism, introduce in new population fittest of old
#   population (if greater than fitness of the individual it is
#   to replace)
    if ielite == 1 && ff(n, newph[1,1]) < fitns[ifit[n]]
        for k = 1 : n
            newph[k, 1] = oldph[k, ifit[np]]
        end
        nnew = nnew - 1
    end

#   replace population
    for i = 1 : np
        for l = 1 : n
            oldph[k, i] = newph[k, i]
        end
#       get fitness using caller's fitness function
        fitns[i]=ff(n,oldph[1,i])        
    end
#   compute new population fitness rank order
    (ifit, jfit) = rnkpop(np,fitns,ifit,jfit)
    return (newph, fitns, ifit, jfit)
end


function encode!(n:: Int, nd:: Int, ph:: Array{Float64, 1}, gn:: Array{Int, 1})
# encode phenotype parameters into integer genotype
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
    
    z = 10.^nd
    ii = 0
    for i = 1 : n
        ip = int(ph[i]*z)
        for j = reverse([1 : nd])
            gn[ii+j] = mod(ip, 10)
        end
        ii = ii + nd
    end
    gn
end

function decode(n:: Int, nd:: Int, gn:: Array{Int, 1})
# decode genotype into phenotype parameters
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
    
    ph = Float64[]

#   z = 10.^(-nd)
#   a^(-b) = 1/ ( a^b)
    z = 1 / (10.^nd)
    ii = 0
    for i = 1 : n
        ip = 0
        for j = 1 : nd
            ip = 10 * ip+gn[ii+j]
        end
#        ph[i]=ip*z
        push!(ph, ip*z)
        ii=ii+nd
    end
    ph
end


function cross!(n:: Int, nd:: Int, pcross:: Float64, gn1:: Array{Int,1}, gn2:: Array{Int, 1})
# breeds two parent chromosomes into two offspring chromosomes
# breeding occurs through crossover. If the crossover probability
# test yields true (crossover taking place), either one-point or
# two-point crossover is used, with equal probabilities.
#
# Compatibility with version 1.0: To enforce 100% use of one-point
# crossover, un-comment appropriate line in source code below

#   Use crossover probability to decide whether a crossover occurs
    if urand() < pcross
#       Compute first crossover point
        ispl = int(urand()*n*nd) + 1
#       Now choose between one-point and two-point crossover 
        if urand() < 0.5
            ispl2 = n*nd
        else
            ispl2 = int(urand()*n*nd) + 1
#           Un-comment following line to enforce one-point crossover
#           ispl2=n*nd
            if ispl2 < ispl
                itmp = ispl2
                ispl2 = ispl
                ispl = itmp
            end
        end
#       Swap genes from ispl to ispl2
        for i = [ispl : ispl2]
            t = gn2[i]
            gn2[i] = gn1[i]
            gn1[i] = t
        end
    end
    gn1, gn2
end

function mutate!(n:: Int, nd:: Int, pmut:: Float64, gn:: Array{Int, 1}, imut:: Int)
# Mutations occur at rate pmut at all gene loci
#   imut=1    Uniform mutation, constant rate
#   imut=2    Uniform mutation, variable rate based on fitness
#   imut=3    Uniform mutation, variable rate based on distance
#   imut=4    Uniform or creep mutation, constant rate
#   imut=5    Uniform or creep mutation, variable rate based on fitness
#   imut=6    Uniform or creep mutation, variable rate based on distance

    fix_it_up = false

#   Decide which type of mutation is to occur
    if imut > 4 && urand() < 0.5
#       CREEP MUTATION OPERATOR
#       Subject each locus to random +/- 1 increment at the rate pmut
        for i = 1 : n
            for j = 1 : nd
                if urand() < pmut
#                   Construct integer
                    loc = (i - 1)*nd + j
                    inc = round( urand() ) * 2 - 1
                    ist = (i - 1)*nd + 1
                    gn[loc] = gn[loc]+inc

                    fix_it_up == false
#                   @printf("%6i %6i %6i", ist, loc, inc)
#                   This is where we carry over the one (up to two digits)
#                   first take care of decrement below 0 case
                    if inc < 0 && gn[loc] < 0
                        if j == 1
                            gn[loc] = 0
                        else
                            for k = reverse([(ist+1) : loc])
                                gn[k] = 9
                                gn[k - 1] = gn[k - 1] - 1
                                if gn[k - 1] > 0
                                    fix_it_up = true
                                    break
                                end
                            end
#                           we popped under 0.00000 lower bound; fix it up
                            if fix_it_up == false
                                if gn[ist] < 0
                                    for l = ist : loc
                                        gn[l] = 0
                                    end
                                end
                            end
                        end
                    end
                    fix_it_up == false
                    if inc > 0 && gn[loc] > 9
                        if j == 1
                            gn[loc] = 9
                        else
                            for k = reverse([(ist+1) : loc])
                                gn[k] = 9
                                gn[k - 1] = gn[k - 1] + 1
                                if gn[k - 1] < 9
                                    fix_it_up = true
                                    break
                                end
                            end
#                           we popped over 9.99999 upper bound; fix it up
                            if fix_it_up == false
                                if gn[ist] > 9
                                    for l = ist : loc
                                        gn[l] = 9
                                    end
                                end
                            end
                        end
                    end
                end
            end # j
        end # i
    else
#   UNIFORM MUTATION OPERATOR
#   Subject each locus to random mutation at the rate pmut
        for i = 1 : (n*nd)
            if urand() < pmut
                gn[i] = int(urand() * 10.0)
            end
        end
    end
    gn
end

function pikaia(ff:: Function, n:: Int, ctrl:: Array{Float64, 1})
# Optimization (maximization) of user-supplied "fitness" function
# ff  over n-dimensional parameter space  x  using a basic genetic
# algorithm method.

# Version 0.0.1   [ 2014 February 21 ]

const NMAX = 52   # NMAX is the maximum number of adjustable parameters (n <= NMAX)
const PMAX = 928  # PMAX is the maximum population (ctrl(1) <= PMAX)
const DMAX = 32   # (default 6) DMAX is the maximum number of Genes (digits) per Chromosome segement (parameter) (ctrl(3) <= DMAX)

#   Init output:
    
    x = rand(n)
    f = 0.
    status  = 0

#   Set control variables from input and defaults
    (status, np, ngen, nd, imut, irep, ielite, ivrb, 
        pcross, pmutmn, pmutmx, pmut, fdif) = setctl(ctrl, n)
    if status != 0
        println(" Control vector (ctrl) argument(s) invalid")
        return (x, f, status)
    end
#   Make sure locally-dimensioned arrays are big enough
    if  n > NMAX || np > PMAX || nd > DMAX 
        println(" Number of parameters, population, or genes too large")
        status = -1
        return (x, f, status)
    end

#  Local variables

    gn1 = [1 : n*nd]
    gn2 = [1 : n*nd]
    fitns = rand(np) 
    ifit = [1 : np]
    jfit = [1 : np]

    ph = rand(n, 2)
    newph = rand(NMAX, PMAX)

    (old_ph, fitns, ifit, jfit) = init_pop(ff, n, np)

#    oldph = Array(Float64, n, np)

#   Compute initial (random but bounded) phenotypes
#    for ip = 1 : np
#        for k = 1 : n
#            oldph[k, ip] = urand()
#        end
#        fitns[ip] = ff(n, oldph[1,ip])
#    end
#   Rank initial population by fitness order
#    (ifit, jfit) = rnkpop(np, fitns)

#   Main Generation Loop
    for ig = 1 : ngen

#       Main Population Loop
        newtot = 0  
        for ip = 1 : np / 2

#           1. pick two parents
            ip1 = select(np,jfit,fdif)
            while true
                ip2 = select(np,jfit,fdif)
                if ip1 != ip2
                    break
                end
            end
#           2. encode parent phenotypes
            encode!(n, nd, oldph[:, ip1], gn1)
            encode!(n, nd, oldph[:, ip2], gn2)

#           3. breed
            cross!(n, nd, pcross, gn1, gn2)
            mutate!(n, nd, pmut, gn1, imut)
            mutate!(n, nd, pmut, gn2, imut)

#           4. decode offspring genotypes 
            ph[:,1] = decode(n, nd, gn1)
            ph[:,2] = decode(n, nd, gn2)
#           5. insert into population
            if irep == 1
                genrep()
            else
                stdrep()
                newtot=newtot+new
            end
        end # End of Main Population Loop

#       if running full generational replacement: swap populations
        if irep = 1
            (old_ph, fitns, ifit, jfit) = new_pop(ff, ielite, ndim, n, np, oldph, fitns)
        end

#       adjust mutation rate?
        if imut == 2 || imut == 3 || imut == 5 || imut == 6
            adjmut()
        end

        if ivrb > 0
            report()
        end

    end # End of Main Generation Loop 

#   Return best phenotype and its fitness
    for k = 1 : n
        x[k] = oldph[k, ifits[np]]
    end
    f = fitns[ifit[np]]
# Return in typle
return (x, f, status)
end # pikaia

end # Pikaia
#
# #3F3434 #DFE8E2
# http://habrahabr.ru/post/161009/
# https://github.com/johnmyleswhite/HopfieldNets.jl
# http://habrahabr.ru/post/125999/
# http://julialang.org/gsoc/2014/

