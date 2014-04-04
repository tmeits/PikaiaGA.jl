# *********************************************************************
# The Julia Language: a fresh approach to technical computing.
# Pikaia is a general purpose function optimization Julia lang module 
# based on a genetic algorithm.
# **********************************************************************
module Pikaia

importall Base

export 
    pikaia, 
    rqsort,
    init_pop,
    new_pop!,
    urand,
    encode!,
    cross!,
    mutate!,
    report,
    adjustment!,
    steady_state_reproduction!,
    select

global _bestft = 0.0
global _pmutpv = 0.0

# *********************************************************************    
function rqsort(n:: Int, a:: Vector{Float64})
# =====================================================================    
# Return integer array p which indexes array a in increasing order.
# Array a is not disturbed.
# ===================================================================== 
    sortperm(a) # the permutation to sort an array
end

# *********************************************************************
function urand()
# =====================================================================    
# Return the next pseudo-random deviate from a sequence which is
# uniformly distributed in the interval [0,1]
# =====================================================================
    
    return rand()
end

function set_ctrl_default()
# Initialize some empty vectors
    ctrl = Float64[]

# First, initialize the random-number generator
    srand(123456)

# Set control variables (use defaults)
    for i= [1:12]
        # Push default values into an array
        ctrl = push!(ctrl, -1.)
    end
    ctrl[1]=1
    ctrl[2]=50.
    
    ctrl
end    


# *********************************************************************
function setctl(ctrl::Array{Float64, 1}, n:: Int)
# =====================================================================    
# Set control variables and flags from input and defaults
# =====================================================================

    const DEFAULT = 
        [100, 500, 5, .85, 2, .005, .0005, .25, 1, 1, 1, 0]

    global _bestft::Float64 = 0.0
    global _pmutpv::Float64 = 0.0

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
#   Print a header
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

# *********************************************************************
function report(ivrb::Int, ndim::Int, n::Int, np::Int, nd::Int, 
    oldph::Matrix{Float64}, fitns::Vector{Float64}, ifit::Vector{Int}, 
    pmut::Float64, ig::Int, nnew::Int)
# =====================================================================    
# Write generation report to standard output
# =====================================================================
    
    global _bestft
    global _pmutpv

    rpt = false

    if pmut != _pmutpv
        _pmutpv = pmut
        rpt = true
    end

    if fitns[ifit[np]] != _bestft
        _bestft = fitns[ifit[np]]
        rpt = true
    end

    if rpt == true || ivrb > 2

#       nint() rounds a real to the nearest result integer
#       Power of 10 to make integer genotypes for display
        ndpwr = iround(10.^nd)

        @printf(" %6i %6i %10.6f %10.6f %10.6f %10.6f \n",
            ig, nnew, pmut, fitns[ifit[np]], fitns[ifit[np-1]], fitns[ifit[np/2]])
        for k=1:n
            @printf(" %10i  %10i  %10i",
                iround(ndpwr*oldph[k,ifit[np]]),
                iround(ndpwr*oldph[k,ifit[np-1]]),
                iround(ndpwr*oldph[k,ifit[np/2]]))
        end
    end
    return true
end    

# call genrep(NMAX,n,np,ip,ph,newph)

# *********************************************************************
function generational_replacement(n::Int, np::Int, ip::Int,
    phenotype:: Vector{Float64})
# =====================================================================    
# full generational replacement: accumulate offspring into new population array
# =====================================================================

    new_phenotype = Float64[]

#   Insert one offspring pair into new population
    i1=2*ip-1
    i2=i1+1
    for k=1:n
        push!(new_phenotype[k, i1], phenotypes[k, 1])
        push!(new_phenotype[k, i2], phenotypes[k, 2])
    end
    new_phenotype
end

# *********************************************************************
function steady_state_reproduction!(ff::Function, ndim::Int, n::Int, 
    np::Int, irep::Int, ielite::Int, ph::Matrix{Float64}, 
    oldph:: Matrix{Float64}, fitns:: Vector{Float64}, 
    ifit::Vector{Int},jfit::Vector{Int})
# =====================================================================    
# steady-state reproduction: insert offspring pair into population
# only if they are fit enough (replace-random if irep=2 or replace-worst if irep=3).
# =====================================================================

    nnew = 0
    goto_j = false
    for j=1:2
#       1. compute offspring fitness (with caller's fitness function)
        fit = ff(n, ph[1, j])
        
#       2. if fit enough, insert in population
        for i=reverse(1:np)
            if fit > fitns[ifit[i]]
#               make sure the phenotype is not already in the population
                if i < np
                    for k=1:n
                        if oldph[k, ifit[i+1]] != ph[k,j]
                            goto_j=true
                            break
                        end
                    end # k
                end # i
                if goto_j == true
                    break
                else
#                   offspring is fit enough for insertion, and is unique
 
#                   (i) insert phenotype at appropriate place in population
                    if irep == 3
                        i1=1
                    elseif ielite == 0 || i == np
                        i1=int(urand()*np)+1
                    else
                        i1=int(urand()*(np-1))+1
                    end
                    if1=ifit[i1]
                    fitns[if1]=fit
                    for k=1:n
                        oldph[k,if1]=ph[k,j]
                    end
#                   (ii) shift and update ranking arrays
                    if i < i1
#                       shift up 
                        jfit[if1]=np-i
                        for k=reverse((i+1):(i1-1))
                            jfit[ifit[k]] = jfit[ifit[k]]-1
                            ifit[k+1]=ifit[k]
                        end
                        ifit[i]=if1
                    else
#                       shift down
                        jfit[if1]=np-i+1
                        for k=(i1+1):i
                            jfit[ifit[k]] = jfit[ifit[k]]+1
                            ifit[k-1]=ifit[k]
                        end
                        ifit[i]=if1
                    end
                    nnew = nnew+1
                    break
                end # goto_j
            end   
        end # i
    end # j
end

# *********************************************************************
function select(population_size::Int, jfit::Vector{Int}, fdif::Float64) # relative fitness differential
# =====================================================================    
# Selects a parent from the population, using roulette wheel
# algorithm with the relative fitnesses of the phenotypes as
# the "hit" probabilities [see Davis 1991, chap. 1].
# =====================================================================   

    idad = 0

    ps1   = population_size+1 
    dice  = urand()*population_size*ps1
    rtfit = 0.0

    for i= 1:population_size
        rtfit = rtfit+ps1+fdif*(ps1-2*jfit[i])
        if rtfit >= dice
            idad=i
            break
        end
    end

    idad

end # select

# *********************************************************************
function  rank_pop(n:: Int, fitnes:: Vector{Float64})
# =====================================================================    
# Rank population by fitness order
# Produce key index and rank order of input array arrin
# =====================================================================

#   Compute the key index
    rank = [1:n] 
    indx = rqsort(n, fitnes)

#   ...and the rank order
    for i = 1:n
        rank[indx[i]] = n-i+1
    end

    return (indx, rank)

end # rank_pop


# *********************************************************************
function init_pop(ff:: Function, n:: Int, np:: Int)
# =====================================================================    
#   Compute initial (random but bounded) phenotypes
# =====================================================================

    old_ph = rand(n, np)
    fitns  = rand(np)

    for ip = 1:np
        for k = 1:n
            old_ph[k, ip] = urand() 
        end
        fitns[ip] = ff(old_ph[:, ip])
    end
#   Rank initial population by fitness order
    (ifit, jfit) = rank_pop(np, fitns)

    return (old_ph, fitns, ifit, jfit)

end # init_pop

# **********************************************************************
function new_pop!(ff::Function, ielite::Int, ndim::Int, n::Int, np::Int, 
    oldph::Matrix{Float64}, newph::Matrix{Float64})
# ======================================================================    
# replaces old population by new; recomputes fitnesses & ranks
# ======================================================================
    
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
    (ifit, jfit) = rank_pop(np,fitns)

    return (newph, fitns, ifit, jfit)

end

# http://whitedwarf.org/metcalfe/node8.htmi
# http://download.hao.ucar.edu/archive/pikaia/

#julia> typeof(10.^6)
#Int64

#julia> typeof(10.0^6)
#Float64
# ??????


# ********************************************************************
function encode!(n::Int, nd::Int, ph::Vector{Float64}, gn::Vector{Int})
# ====================================================================    
# encode phenotype parameters into integer genotype
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
# ====================================================================

    z  = 10.0^nd
    ii = 0
    print("z= "); println(z)

    for i = 1:n
        ip = int(ph[i]*z)
        @printf("ip = %6i\n", ip)
        for j = reverse([1:nd])
#            print(ii+j)
#            println(mod(ip, 10))
            gn[ii+j] = mod(ip, 10)
            ip=int(floor(ip/10))
            @printf("%6i %6i %9.4f %6i \n", ii, j,  mod(ip, 10), ip)
        end
        ii = ii+nd
    end

    gn
end

# *********************************************************************
function decode(n::Int, nd::Int, gn::Vector{Int})
# =====================================================================    
# decode genotype into phenotype parameters
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
# =====================================================================
    
    ph = Float64[]

#   z = 10.^(-nd)
#   a^(-b) = 1/ ( a^b)
    z  = 1./(10.^nd)
    ii = 0

    for i = 1:n
        ip = 0
        for j = 1:nd
            ip = 10*ip+gn[ii+j]
        end
#        ph[i]=ip*z
        push!(ph, ip*z)
        ii=ii+nd
    end

    ph
end

# *********************************************************************
function cross!(n::Int, nd::Int, pcross::Float64, gn1::Vector{Int},
    gn2::Vector{Int})
# =====================================================================    
# breeds two parent chromosomes into two offspring chromosomes
# breeding occurs through crossover. If the crossover probability
# test yields true (crossover taking place), either one-point or
# two-point crossover is used, with equal probabilities.
#
# Compatibility with version 1.0: To enforce 100% use of one-point
# crossover, un-comment appropriate line in source code below
# =====================================================================

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

# *********************************************************************
function mutate!(n::Int, nd::Int, pmut::Float64, gn::Vector{Int}, imut::Int)
# =====================================================================    
# Mutations occur at rate pmut at all gene loci
#   imut=1    Uniform mutation, constant rate
#   imut=2    Uniform mutation, variable rate based on fitness
#   imut=3    Uniform mutation, variable rate based on distance
#   imut=4    Uniform or creep mutation, constant rate
#   imut=5    Uniform or creep mutation, variable rate based on fitness
#   imut=6    Uniform or creep mutation, variable rate based on distance
# =====================================================================
    fix_it_up = false

#   Decide which type of mutation is to occur
    if imut > 4 && urand() < 0.5
#       CREEP MUTATION OPERATOR
#       Subject each locus to random +/- 1 increment at the rate pmut
        for i = 1:n
            for j = 1:nd
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
                            for k = reverse([(ist+1):loc])
                                gn[k] = 9
                                gn[k-1] = gn[k-1]-1
                                if gn[k-1] > 0
                                    fix_it_up = true
                                    break
                                end
                            end
#                           we popped under 0.00000 lower bound; fix it up
                            if fix_it_up == false
                                if gn[ist] < 0
                                    for l = ist:loc
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
                            for k = reverse([(ist+1):loc])
                                gn[k] = 9
                                gn[k-1] = gn[k-1]+1
                                if gn[k-1] < 9
                                    fix_it_up = true
                                    break
                                end
                            end
#                           we popped over 9.99999 upper bound; fix it up
                            if fix_it_up == false
                                if gn[ist] > 9
                                    for l = ist:loc
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
        for i = 1:(n*nd)
            if urand() < pmut
                gn[i] = int(urand() * 10.0)
            end
        end
    end
    gn
end

# *********************************************************************
function adjustment!(ndim::Int, n::Int, np::Int, oldph::Matrix{Float64}, 
    fitns::Vector{Float64}, ifit::Vector{Int}, pmutmn::Float64, 
    pmutmx::Float64, pmut::Float64, imut::Int) 
# =====================================================================    
# dynamical adjustment of mutation rate;
#    imut=2 or imut=5 : adjustment based on fitness differential
#                       between best and median individuals
#    imut=3 or imut=6 : adjustment based on metric distance
#                       between best and median individuals

#     dynamical adjustment of mutation rate;
#        imut=2 or imut=5 : adjustment based on fitness differential
#                           between best and median individuals
#        imut=3 or imut=6 : adjustment based on metric distance
#                           between best and median individuals
# =====================================================================

    const rdiflo=0.05
    const rdifhi=0.25
    const delta=1.5    
    
    if imut == 2 || imut == 5
#   Adjustment based on fitness differential
        rdif=abs(fitns[ifit[np]]-fitns[ifit[np/2]])/
                (fitns[ifit[np]]+fitns[ifit[np/2]])
    elseif imut == 3 || imut == 6 
#   Adjustment based on normalized metric distance        
        rdif=0.0
        for i=1:n
            rdif=rdif+( oldph[i,ifit[np]]-oldph[i,ifit[np/2]] )^2
        end
        rdif=sqrt(rdif)/float(n)
    end

    if rdif < rdiflo 
        pmut=min(pmutmx, pmut*delta)
    elseif rdif > rdifhi
        pmut=max(pmutmn, pmut/delta)
    end
end

# *********************************************************************
function pikaia(ff::Function, n::Int, ctrl::Vector{Float64})
# =====================================================================
# Optimization (maximization) of user-supplied "fitness" function
# ff  over n-dimensional parameter space  x  using a basic genetic 
# algorithm method.
# =====================================================================

# Version 0.0.1   [ 2014 February 21 ]

#   Output: 
    x = rand(n); f = 0.0; status  = 0
#
#    o Array  x[1:n]  is the "fittest" (optimal) solution found,
#       i.e., the solution which maximizes fitness function ff
#
#    o Scalar  f  is the value of the fitness function at x
#
#    o Integer  status  is an indicator of the success or failure
#       of the call to pikaia (0=success; non-zero=failure)
#
#

#   Constants
    const NMAX = 32  
    const PMAX = 328  
    const DMAX = 6  
#
#    o NMAX is the maximum number of adjustable parameters
#      (n <= NMAX)
#
#    o PMAX is the maximum population (ctrl(1) <= PMAX)
#
#    o DMAX is the maximum number of Genes (digits) per Chromosome
#      segement (parameter) (ctrl(3) <= DMAX)
#

#   Set control variables from input and defaults
    (status, np, ngen, nd, imut, irep, ielite, ivrb, 
        pcross, pmutmn, pmutmx, pmut, fdif) = setctl(ctrl, n)
    if status != 0
        warn(" Control vector (ctrl) argument(s) invalid")
        return (x, f, status)
    end

#   Local variables
    gn1   = [1:n*nd]
    gn2   = [1:n*nd]
    fitns = rand(np) 
    ifit  = [1:np]
    jfit  = [1:np]

    ph    = rand(n, 2)
    newph = rand(n, np)
 
#   Make sure locally-dimensioned arrays are big enough
    if  n > NMAX || np > PMAX || nd > DMAX 
        warn(" Number of parameters, population, or genes too large")
        status = -1
        return (x, f, status)
    end

#   Compute initial (random but bounded) phenotypes
#   Rank initial population by fitness order

    (old_ph, fitns, ifit, jfit) = init_pop(ff, n, np)

   
    oldph = Array(Float64, n, np)

#   Main Generation Loop
    for ig = 1:ngen

#       Main Population Loop
        newtot = 0  
        for ip = 1:(np/2)

#           1. pick two parents
            ip1 = select(np, jfit, fdif)
            while true
                ip2 = select(np, jfit, fdif)
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
        if irep == 1
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
# http://omega.sp.susu.ac.ru/books/conference/PaVT2008/papers/Short_Papers/030.pdf
# http://www.hao.ucar.edu/modeling/pikaia/
# http://getbootstrap.com/examples/carousel/#
