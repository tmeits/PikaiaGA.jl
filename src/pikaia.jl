# *********************************************************************
# Pikaia is a general purpose function optimization Julia lang module 
# based on a genetic algorithm.
# **********************************************************************

# include("pikaia.jl")
# include("testfunctions.jl")

module Pikaia

importall Base

export

# SYSTEM FUNCTION
    rqsort,
    init_phenotype,
    init_pop,
    new_pop!,
    urand,
    report,
    steady_state_reproduction!,
    select,
    select2
    
# GENETICS MODULE

    encode,
    decode,
    cross,
    one_point_crossover,
    mutate,
    adjust_mutation,

# MAIN FUNCTION    
  
    pikaia

global _bestft = 0.0
global _pmutpv = 0.0
global _tprint = false

function dbg(str::String)
    global _tprint

    if _tprint == true
        print(str)
    end

    true
end    

# *********************************************************************
function fortran_int(var:: Float64)
# =====================================================================
# integer part var:: Float64
# =====================================================================    

    int_var = int(floor(var))

    return int_var
end    

# *********************************************************************    
function rqsort(n:: Int, a:: Vector{Float64})
# =====================================================================    
# Return integer array p which indexes array a in increasing order.
# Array a is not disturbed.
# ===================================================================== 

    return sortperm(a) # the permutation to sort an array
end

# *********************************************************************
function urand()
# =====================================================================    
# Return the next pseudo-random deviate from a sequence which is
# uniformly distributed in the interval [0,1]
# =====================================================================
    
    return rand()
end

# *********************************************************************
function set_ctrl_default(seed:: Int)
# =====================================================================    
# Initialize some empty vectors
# =====================================================================

    const size_ctrl = 12

    index = [1: size_ctrl]
    ctrl  = Float64[]

# First, initialize the random-number generator
    srand(seed)

# Set control variables (use defaults)
    for i in index
        # Push default values into an array
        ctrl = push!(ctrl, -1.0)
    end

    ctrl[1] = 1
    ctrl[2] = 50.0
    
    return ctrl
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

    for i = 1:length(ctrl)
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
function report(
    ivrb::Int,              # 
    ndim::Int,              #
    n::Int,                 #
    np::Int,                #
    nd::Int,                #
    oldph::Matrix{Float64}, #
    fitns::Vector{Float64}, #
    ifit::Vector{Int},      #
    pmut::Float64,          #
    ig::Int,                #
    nnew::Int               #
    )
# =====================================================================    
# Write generation report to standard output
# =====================================================================
    
    global _bestft
    global _pmutpv

    rpt = false

    if pmut != _pmutpv
        _pmutpv = pmut
        rpt     = true
    end

    if fitns[ifit[np]] != _bestft
        _bestft = fitns[ifit[np]]
        rpt     = true
    end

    if rpt == true || ivrb > 2

#       nint() rounds a real to the nearest result integer
#       Power of 10 to make integer genotypes for display
        ndpwr = iround(10.^nd)

        @printf(" %6i %6i %10.6f %10.6f %10.6f %10.6f \n",
            ig, nnew, pmut, fitns[ifit[np]], fitns[ifit[np-1]], fitns[ifit[np/2]])
        for k = 1:n
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
function generational_replacement(
    n::Int,  # n is the number of adjustable parameters
    np::Int, # number of individuals in a population
    ip::Int, # Population Loop
    ph:: Array{Float64,2}, # 
    new_phenotype:: Matrix{Float64} #
    ) 
# =====================================================================    
# full generational replacement: accumulate offspring into new population array
# Inserts offspring into population, for full generational replacement
# =====================================================================

    local_new_phenotype = copy(new_phenotype)

#   Insert one offspring pair into new population
    i1 = 2*ip-1
    i2 = i1+1

    for k=1:n
        local_new_phenotype[k, i1] = ph[k, 1]
        local_new_phenotype[k, i2] = ph[k, 2]
    end
    
    return local_new_phenotype
end 
# http://qai.narod.ru/GA/meta-heuristics_3.pdf page 17

#call stdrep(ff,NMAX,n,np,irep,ielite,
#     +                     ph,oldph,fitns,ifit,jfit,new)
# *********************************************************************
function steady_state_reproduction!(
    ff::Function,           # fitness function
#   ndim::Int,              # (constant) MAX is the number of adjustable parameters
    n::Int,                 # n is the number of adjustable parameters
    np::Int,                # number of individuals in a population
    irep::Int,              # ctrl(10) type replacement 1,2,3
    ielite::Int,            # ctrl(11) elitism flag; 0/1=off/on (default is 0)
                            # (Applies only to reproduction plans 1 and 2)

    ph::Matrix{Float64},    #
    oldph:: Matrix{Float64},#
    fitns:: Vector{Float64},#
    ifit::Vector{Int},      #
    jfit::Vector{Int}       #
    )
# =====================================================================    
# steady-state reproduction: insert offspring pair into population
# only if they are fit enough (replace-random if irep=2 or replace-worst if irep=3).
# =====================================================================
# http://www.epathtrade.com/product.asp?/G-5600E-1D-G-SHOCK-Tough-Solar-Sport-Watch,194.htm
# I saw in a network the management for Julia of the Google's R Style Guide type, and now I can't find. Throw the link.
    nnew = 0
    goto_j = false

    for j = 1:2

#       1. compute offspring fitness (with caller's fitness function)
        fit = ff(ph[:, j])
        
#       2. if fit enough, insert in population
        for i = reverse(1:np)
            if fit > fitns[ifit[i]]
                @printf("steady_state_reproduction! fit>")
#               make sure the phenotype is not already in the population
                if i < np
                    for k = 1:n
                        if oldph[k, ifit[i+1]] != ph[k,j]
                            goto_j = true
                            break
                        end
                    end # k
                end # if
                if goto_j == true
#                   offspring is fit enough for insertion, and is unique
 
#                   (i) insert phenotype at appropriate place in population
                    if irep == 3
                        i1 = 1
                    elseif ielite == 0 || i == np
                        i1 = fortran_int((urand()*np))+1
                    else
                        i1 = fortran_int((urand()*(np-1)))+1
                    end

                    if1        = ifit[i1]
                    fitns[if1] = fit
                    
                    for k = 1:n
                        oldph[k,if1] = ph[k,j] # 
                    end

#                   (ii) shift and update ranking arrays
                    if i < i1
#                       shift up 
                        jfit[if1] = np-i

                        for k = reverse((i+1):(i1-1))
                            jfit[ifit[k]] = jfit[ifit[k]]-1
                            ifit[k+1]     = ifit[k]
                        end
                        
                        ifit[i] = if1
                    else
#                       shift down
                        jfit[if1] = np-i+1

                        for k = (i1+1):i
                            jfit[ifit[k]] = jfit[ifit[k]]+1
                            ifit[k-1]     = ifit[k]
                        end

                        ifit[i] = if1
                    end

                    nnew = nnew+1
                    break

                end # if goto_j
            end  # if fit 
        end # for i
    end # for j
    
    return oldph

end # function steady_state_reproduction!
# http://blog.yhathq.com/posts/julia-neural-networks.html
# http://arstechnica.com/science/2014/05/scientific-computings-future-can-any-coding-language-top-a-1950s-behemoth/
# http://juliawebstack.org/

# *********************************************************************
function select(population_size::Int, jfit::Vector{Int}, fdif::Float64) 
# relative fitness differential
# =====================================================================    
# Selects a parent from the population, using roulette wheel
# algorithm with the relative fitnesses of the phenotypes as
# the "hit" probabilities [see Davis 1991, chap. 1].
# =====================================================================   

    idad  = 0

    ps1   = population_size+1 
    dice  = urand()*population_size*ps1        # faites vos jeux...
    rtfit = 0.0

    for i = 1:population_size                  # this is the running sum
        rtfit = rtfit+ps1+fdif*(ps1-2*jfit[i]) # compute relative
        if rtfit >= dice                       # fitness on the fly
            idad = i                           # a parent has been found
            break
        end
    end

    return idad
end #select

# for i=1:10 if i == 5 println(i); break end end

# *********************************************************************
function select2(population_size::Int, jfit::Vector{Int}, fdif::Float64)

    ip1 = -1
    ip2 = -1

#   1. pick two parents
    ip1 = select(population_size, jfit, fdif)
#    @printf("selct2::ip1= %6i\n", ip1)
    while true
        ip2 = select(population_size, jfit, fdif)
#        @printf("selct2::ip2= %6i\n", ip2)
        if ip1 != ip2
            break
        end
    end
    
    return (ip1, ip2)
end

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
function init_phenotype(number_parameters:: Int, population_size:: Int)
# =====================================================================
# Compute initial (random but bounded) phenotypes
# =====================================================================
    
    phenotypes = Array(Float64, number_parameters, population_size)

    for i = 1:population_size
        for j = 1:number_parameters
            phenotypes[j, i] = rand()
        end
    end

    return phenotypes
end
# *********************************************************************
function init_pop(ff:: Function, n:: Int, np:: Int)
# =====================================================================    
# Compute initial (random but bounded) phenotypes
# =====================================================================

    old_ph = Array(Float64, n, np)
    fitns  = Array(Float64, np)

    for ip = 1:np                      # np individuals in population
        for k = 1:n                    # n parameters per individual
            old_ph[k, ip] = urand()    # Random initialization
        end
        fitns[ip] = ff(old_ph[:, ip])  # Compute fitness
    end
#   Rank initial population by fitness order
    (ifit, jfit) = rank_pop(np, fitns) #

    return (old_ph, fitns, ifit, jfit)

end #init_pop

# **********************************************************************
function new_pop!(
    ff::Function,           #
    ielite::Int,            #
    n::Int,                 #
    np::Int,                #

    oldph::Matrix{Float64}, #
    newph::Matrix{Float64}  #
    )
# ======================================================================    
# replaces old population by new; recomputes fitnesses & ranks
# ======================================================================
    
    nnew  = np
    newph = copy(oldph)

#   if using elitism, introduce in new population fittest of old
#   population (if greater than fitness of the individual it is to replace)
    if ielite == 1 && ff(n, newph[1,1]) < fitns[ifit[n]]
        for k = 1:n
            newph[k, :] = oldph[k, ifit[np]]
        end
        nnew = nnew - 1
    end

#   replace population
    for i = 1:np
        for l = 1:n
            oldph[k, i] = newph[k, i]
        end
#       get fitness using caller's fitness function
        fitns[i] = ff(n, oldph[: , i])        
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

#[16:55] <tmeits> hi
#[16:56] <tmeits> julia> typeof(10.^6) Int64  julia> typeof(10.0^6) Float64 ??????
#[17:01] <acetoline> 10.^6 is parsed as 10 .^ 6
#[17:01] <acetoline> .^ is an operator
#[17:02] <acetoline> you probably mean 10. ^ 6
#[17:02] <tmeits> thanks is clear
#[17:03] <tmeits> I tried precisely once again thanks

# using ASCIIPlots

# ********************************************************************
function encode(n::Int, nd::Int, ph::Vector{Float64})
# ====================================================================    
# encode phenotype parameters into integer genotype
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
# ====================================================================

    gn = Array(Int, nd)

    z  = 10.0 ^ nd
    ii = 0

    for i = 1:n
        ip = int(ph[i]*z)
#       @printf("ip = %6i\n", ip)
        for j = reverse([1:nd])
            gn[ii+j] = mod(ip, 10)
            ip = fortran_int((ip/10))
#           @printf("%6i %6i %9.4f %6i \n", ii, j,  mod(ip, 10), ip)
        end
        ii = ii+nd
    end

    return gn
end

# *********************************************************************
function decode(n::Int, nd::Int, gn::Vector{Int})
# =====================================================================    
# decode genotype into phenotype parameters
# ph(k) are x,y coordinates [ 0 < x,y < 1 ]
# =====================================================================
    
    ph = Float64[]    
#    ph = Array(Float64, n)

#   z = 10.^(-nd)
#   a^(-b) = 1/ ( a^b)
    z  = 1.0/(10.0 ^ nd)
    ii = 0

    for i = 1:n
        ip = 0
        for j = 1:nd
            ip = 10*ip+gn[ii+j]
        end
#        ph[i]=ip*z
        push!(ph, ip*z)
        ii = ii+nd
    end

    return ph
end

# *********************************************************************
function get_random_int(rand_num_min::Int, rand_num_max::Int)
# =====================================================================    
# Returns a random integer between min and max
# =====================================================================

    rand_int = (fortran_int((rand() * (rand_num_max - rand_num_min + 1)) + rand_num_min))

    rand_int
end

# *********************************************************************
function get_random_int(num::Int, rand_num_min::Int, rand_num_max::Int)
# =====================================================================    
# Returns a vector casual whole within a minimum and a maximum
# =====================================================================

    rand_int = Int[]

    for i=1:num
        push!(rand_int, 
            (fortran_int((rand() * (rand_num_max - rand_num_min + 1)) + rand_num_min)))
    end

    rand_int
end


# http://mathmod.aspu.ru/images/File/ebooks/GAfinal.pdf 
# *********************************************************************
function one_point_crossover(n::Int, nd::Int, pcross::Float64,
    gn1::Vector{Int}, gn2::Vector{Int})
# =====================================================================  
# breeds two parent chromosomes into two offspring chromosomes
# breeding occurs through crossover. 
# =====================================================================

    gen1 = Int[]
    gen2 = Int[]

    for i = 1:length(gn1)
        push!(gen1, gn1[i])
        push!(gen2, gn2[i])
    end
#     gen1 = copy(gn1)
#     gen2 = copy(gn2)
    ce = false

    if urand() < pcross
        ce   = true
        ispl = fortran_int((urand()*n*nd))+1 # choose cutting point
#       @printf("%7i\n",ispl)
        for i = ispl:n*nd
            t       = gen2[i]
            gen2[i] = gen1[i]
            gen1[i] = t
        end
    end

    return (ce, gen1, gen2)
end

# *********************************************************************
function cross!(n::Int, nd::Int, pcross::Float64,
    gn1::Vector{Int}, gn2::Vector{Int})
# =====================================================================    
# breeds two parent chromosomes into two offspring chromosomes
# breeding occurs through crossover. If the crossover probability
# test yields true (crossover taking place), either one-point or
# two-point crossover is used, with equal probabilities.
#
# Compatibility with version 1.0: To enforce 100% use of one-point
# crossover, un-comment appropriate line in source code below
# =====================================================================

    ce = false

#   Use crossover probability to decide whether a crossover occurs
    if urand() < pcross
        ce = true
        dbg("cross! yes!\n")
#       Compute first crossover point
        ispl = fortran_int((urand()*n*nd))+1
#       Now choose between one-point and two-point crossover 
        if urand() < 0.5
            ispl2 = n*nd
        else
            ispl2 = fortran_int((urand()*n*nd))+1
#           Un-comment following line to enforce one-point crossover
#           ispl2=n*nd
            if ispl2 < ispl
                itmp  = ispl2
                ispl2 = ispl
                ispl  = itmp
            end
        end
#       Swap genes from ispl to ispl2
        for i = [ispl:ispl2]
            t = gn2[i]
            gn2[i] = gn1[i]
            gn1[i] = t
        end
    end

    return ce, gn1, gn2
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
    me        = false

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
        for i = 1:n*nd
            if urand() < pmut
                gn[i] = fortran_int((urand()*10.0))
            end
        end
    end # if

    gn
end

# *********************************************************************
function adjust_mutation(
    n::Int,                 #
    np::Int,                #
    oldph::Matrix{Float64}, #
    fitns::Vector{Float64}, #
    ifit::Vector{Int},      #
    pmutmn::Float64,        #
    pmutmx::Float64,        #
    pmut::Float64,          #
    imut::Int               #
    ) 
# =====================================================================    
# dynamical adjustment of mutation rate;
#    imut=2 or imut=5 : adjustment based on fitness differential
#                       between best and median individuals
#    imut=3 or imut=6 : adjustment based on metric distance
#                       between best and median individuals

    const rdiflo = 0.05
    const rdifhi = 0.25
    const delta  = 1.5   

    rdif = 0.5
    
    if imut == 2 || imut == 5
#   Adjustment based on fitness differential
        rdif = abs(fitns[ifit[np]]-fitns[ifit[np/2]])/
                  (fitns[ifit[np]]+fitns[ifit[np/2]])
    elseif imut == 3 || imut == 6 
#   Adjustment based on normalized metric distance        
        rdif = 0.0
        for i = 1:n
            rdif = rdif+( oldph[i,ifit[np]]-oldph[i,ifit[np/2]] )^2
        end
        rdif = sqrt(rdif)/float(n)
    end

    if rdif < rdiflo 
        pmut = min(pmutmx, pmut*delta)
    elseif rdif > rdifhi
        pmut = max(pmutmn, pmut/delta)
    end

    return pmut 
end
# http://www.raai.org/resurs/papers/kolomna2009/doklad/Sinyuk_Akopov.doc



# *********************************************************************
function pikaia(ff::Function, n::Int, ctrl::Vector{Float64})
# =====================================================================
# Optimization (maximization) of user-supplied "fitness" function
# ff over n-dimensional parameter space  x  using a basic genetic 
# algorithm method.
# =====================================================================

# Version 0.0.1   [ 2014 February 21 ]

#   Output:
    x       = rand(n); 
    f       = 0.0; 
    status  = 0
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

#   Constants !!!
    const NMAX = 32  
    const PMAX = 328  
    const DMAX = 6  
#
#   o NMAX is the maximum number of adjustable parameters
#      (n <= NMAX)
#
#   o PMAX is the maximum population 
#       (ctrl(1) <= PMAX)
#
#   o DMAX is the maximum number of Genes (digits) per Chromosome
#       segement (parameter) 
#       (ctrl(3) <= DMAX)
#

#   Set control variables from input and defaults
    (status, np, ngen, nd, imut, irep, ielite, ivrb, 
        pcross, pmutmn, pmutmx, pmut, fdif) = setctl(ctrl, n)

    if status != 0
        warn(" Control vector (ctrl) argument(s) invalid")

        return (x, f, status)
    end

#   Local variables
    gn1    = Array(Int, n*nd)
    gn2    = Array(Int, n*nd)
    fitns  = Array(Float64, np) 
    ifit   = Array(Int, np)
    jfit   = Array(Int, np)

    ph     = Array(Float64, n, 2)
    new_ph = Array(Float64, n, np)
 
#   Make sure locally-dimensioned arrays are big enough
    if  n > NMAX || np > PMAX || nd > DMAX 
        warn(" Number of parameters, population, or genes too large")
        status = -1

        return (x, f, status)
    end

#   Compute initial (random but bounded) phenotypes
#   Rank initial population by fitness order

    (old_ph, fitns, ifit, jfit) = init_pop(ff, n, np)
 
#    oldph = Array(Float64, n, np)

#   Main Generation Loop
    for ig = 1:ngen

#       Main Population Loop
        newtot = 0  
        for ip = 1:int(np/2)

#           selection strategies            
#           1. pick two parents
            ip1 = select(np, jfit, fdif)
            while true
                ip2 = select(np, jfit, fdif)
                if ip1 != ip2
                    break
                end
            end
#           2. encode parent phenotypes
            encode!(n, nd, old_ph[:, ip1], gn1)
            encode!(n, nd, old_ph[:, ip2], gn2)

#           3. breed
            cross!(n, nd, pcross, gn1, gn2)
            mutate!(n, nd, pmut, gn1, imut)
            mutate!(n, nd, pmut, gn2, imut)

#           4. decode offspring genotypes 
            ph[:,1] = decode(n, nd, gn1)
            ph[:,2] = decode(n, nd, gn2)
            
#           5. insert into population
# irep == reproduction plan; 1/2/3=Full generationalreplacement/Steady-state-replace-random/
# Steady-state-replace-worst (default is 3)

            if irep == 1
# full generational replacement: accumulate offspring into new population array              
                new_ph = generational_replacement(n, np, ip, ph, new_ph)
            else
                steady_state_reproduction!(ff, n, np, irep, ielite, ph, oldphm, fitns, ifit, jfit)
                newtot = newtot+new
            end # insertion/storage completed

        end # End of Main Population Loop

#       if running full generational replacement: swap populations
        if irep == 1
            (old_ph, fitns, ifit, jfit) = 
                new_pop(ff, ielite, ndim, n, np, old_ph, fitns)
        end

#       adjust mutation rate?
        if imut == 2 || imut == 3 || imut == 5 || imut == 6
            adjmut(NMAX, n, np, oldph, fitns, ifit, pmutmn, pmutmx, pmut, imut)
        end

# printed output 0/1/2=None/Minimal/Verbose (default is 0)

        if ivrb > 0
            report(ivrb, NMAX, n, np, nd, old_ph, fitns, ifit, pmut, ig, newtot)
        end

    end # End of Main Generation Loop 

#   Return best phenotype and its fitness
    for k = 1:n
        x[k] = old_ph[k, ifits[np]]
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
