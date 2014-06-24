# «triple-triple redundant architecture»

using Base.Test
using TestFunctions
using Pikaia

a = rand(1 : 50)

@test sort(a) == a[Pikaia.rqsort(length(a), a, [1:50])]

@test Pikaia.setctl([100, 500, 7, .85, 0, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 5
@test Pikaia.setctl([100, 500, 7, .85, 7, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 5
@test Pikaia.setctl([100, 500, 7, .85, 2, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 0

ng = [1: 55*6]
@test length(Pikaia.encode(55, 6, rand(55), ng) ) == 330
@test Pikaia.encode(55, 6, rand(55), ng)  == ng
ng = [1: 55*6]
@test Pikaia.encode(55, 6, rand(55), ng)  !=  [1: 55*6]

@test Pikaia.mutate!(55,9, 0.85, [1 : 55*9], 1) != [1 : 55*9]

pop_ph=rand(50,850);
gn=Pikaia.encode!(50,6,pop_ph[:,1],[1:50*6]);
Pikaia.decode(50, 5, gn);

# number of output arguments
@test length(TestFunctions.rastriginsfcn([0., 0.])) == 1

# type of output
@test typeof(TestFunctions.rastriginsfcn([0., 0.])) == Float64

# result
@test TestFunctions.rastriginsfcn([0., 0.]) == 0.


@test TestFunctions.rastrigin([0.]) == 0.0
@test TestFunctions.rastrigin([0., 0.]) == 0.0
@test TestFunctions.rastrigin(rand(10)*0.0) == 0.0

 Pikaia.init_pop(TestFunctions.rastrigin,1,4)
 Pikaia.set_ctrl_default()
 Pikaia.select(4, Pikaia.init_pop(TestFunctions.rastrigin,1,4)[4],1.)

(ph,fitns,ifit,jfit)=Pikaia.init_pop(TestFunctions.rastrigin,1,4)
g1=[1:5]
Pikaia.encode!(1,5,ph[:,1],g1)

for ss=1:32 g1=[1:ss]; println(Pikaia.encode!(1,ss,ph[:,3],g1)) ;print("\n")end

for i=1:1000
    val=rand(); 
    nd=int(rand()*32); 
    Pikaia.decode(1,nd,Pikaia.encode!(1,nd,[val],[1:nd]))[1] == val;
end
val=rand(); nd=6; (Pikaia.decode(1,nd,Pikaia.encode!(1,nd,[val],[1:nd]))[1], val)


@test Pikaia.encode(1,6,[0.27749011]) == [2, 7, 7, 4, 9, 0]

function enc_dec(n::Int, d::Int)
    for i =1:n 
        val = rand()
        nd  = d #int(floor(rand()*16))+1
            
        v1  = @sprintf("%12.10f", Pikaia.decode(1, nd,
            Pikaia.encode(1, nd, [val]))[1])
        v2  = @sprintf("%12.10f", val)

        @test v1[1:(d-2)] == v2[1:(d-2)]
    end
    true
end

enc_dec(10000,5)

include("pikaia.jl")
Pikaia.encode!(2,8, [0.123456789, 0.987654321],[1:8*2])
Pikaia.decode(1,8,[1,2,3,4,5,6,7,8])
Pikaia.decode(1,8,[1,2,3,4,5,6,7,9])
Pikaia.decode(1,8,[9,8,7,6,5,4,3,2])

gn1=Pikaia.encode!(2,8, [0.123456789, 0.987654321],[1:8*2])
Pikaia.cross!(2,8,0.15,gn1,reverse(gn1))

r=map((x)->Pikaia.get_random_int(-16,16),[1:10000]);
rt=map((x)-> x >= -16 && x <= 16, r);
@test sim(rt) == 10000

function test_cross()
# include in functions doesn't work! ! !
# at change in the loaded file of change don't come into force    
#   include("pikaia.jl")

    all_true = Bool[]
    res_true = true

    all_true2 = Bool[]
    res_true2 = true

    for i=1:10000
        (istr, a1, a2)=Pikaia.one_point_crossover(1,8,0.67,[1:8],[1:8])
        push!(all_true, a1 == a2)
    end
    for i=1:10000
        if all_true[i] != true
            res_true = false
        end
    end

    for i=1:10000
        gn1=Pikaia.get_random_int(10,0,8)
        gn2=Pikaia.get_random_int(10,0,8)
        (tr, g1, g2) =Pikaia.one_point_crossover(1,10,0.67,gn1,gn2)
        if tr == true
            push!(all_true2, g1 != gn1 && g2 != gn2 )
        else
            push!(all_true2, true)
        end
    end
    for i=1:10000
        if all_true2[i] != true
            res_true2 = false
        end
    end

    res_true, res_true2 == false
end   

function test_cross!()

    const inum=9

    all_true = Bool[]
    res_true = true

    for i=1:10000
        gn1=Pikaia.get_random_int(10,0,inum)
        gn2=Pikaia.get_random_int(10,0,inum)
     
        gc1=copy(gn1); gc2=copy(gn2); 

        (tr, g1, g2) = Pikaia.cross!(1,10,rand(),gn1,gn2)

        if tr == true
            push!(all_true, gc1 != gn1 && gc2 != gn2 )
        end
    end
    for i=1:length(all_true)
        if all_true[i] == false
            res_true = false
        end
    end

    res_true

end   

Pikaia.mutate!(1,10,0.5,Pikaia.get_random_int(10,0,10),6)
i=0
for i=1:100000
    rnd  = Pikaia.get_random_int(10,0,10)
    temp = copy(rnd)
    if Pikaia.mutate!(1,10,rand(),rnd,
        Pikaia.get_random_int(1,6)) != temp
        i=i+1
    end
end
i


include("pikaia.jl")
include("testfunctions.jl")
ph=Pikaia.init_phenotype(1, 10)
(ph,fitns,ifit,jfit)=Pikaia.init_pop(TestFunctions.rastrigin,1,10)
(ip1, ip2)=Pikaia.select2(10, ifit, 0.7)
ph2=Array(Float64,1,2)
ph2[:,1]=ph[:,ip1]
ph2[:,2]=ph[:,ip2]
Pikaia.generational_replacement(1,10,4,ph2, ph)

# http://www.reddit.com/r/Julia/
for idx=1:1000
    tic()
    (ph,fitns,ifit,jfit)=Pikaia.init_pop(TestFunctions.rastrigin,1,10)
    (ip1, ip2)=Pikaia.select2(10, ifit, 0.7)
    ph2=Array(Float64,1,2)
    ph2[:,1]=ph[:,ip1]
    ph2[:,2]=ph[:,ip2]
    Pikaia.steady_state_reproduction!(TestFunctions.rastrigin,1,10,2,0,ph2,ph,fitns,ifit,jfit)
    toc()
end

(ph,fitns,ifit,jfit)=Pikaia.init_pop(TestFunctions.rastrigin,1,10)
Pikaia.adjust_mutation(1,6, ph, fitns, ifit, 0.15, 0.15, 0.15, 2)

# TODO create function steady_state_reproduction_test

# http://www.jstatsoft.org/v53/i04/paper
# Examples
# Function optimization on one dimension


function test_rescaling()

    tr = true

    for i = 1:100000
        min_max = sort(Pikaia.get_random_int(2, -100 ,100)*1.)
        while true
            if  min_max[1] == min_max[2]
                min_max = sort(Pikaia.get_random_int(2, -100 ,100)*1.)
            else
                break
            end
        end
        rnd     = rand()
        res     = Pikaia.rescaling(rnd ,0., 1., min_max[1], min_max[2])
        rnd_    = Pikaia.rescaling(res ,min_max[1], min_max[2], 0., 1.)
        #@printf("min= %9.4f max= %9.4f res= %9.4f rnd = %9.4f == %9.4f\n",
        #    min_max[1], min_max[2], res, rnd, rnd_)
        if strip(@sprintf("%9.4f", rnd)) != strip(@sprintf("%9.4f", rnd_))
            tr = false
        end
    end

    return tr
end

@test test_rescaling()

function ff_easy(x)

    return abs(x[1]) + cos(x[1])
end  

function ff_rescaling(x)
#   appropriately rescaling x
    sx = x[1]
    rx = Pikaia.rescaling(sx,0.,1., -20., 20.) 
#   search min function
    return -1.0*(abs(rx) + cos(rx))
end   

function ff_rescaling(x, rmin::Float64, rmax::Float64)
#   appropriately rescaling x
    rx = Pikaia.rescaling(x,0.,1., rmin, rmax) 

    return abs(rx[1]) + cos(rx[1])
end   

using ASCIIPlots
lineplot([-20:20], map(ff_easy,[-20:20]))
scatterplot([-20:20], map(ff_easy,[-20:20]))

#=
	-------------------------------------------------------------
	|\                                                          /| 20.41
	| \                                                       /  |
	|  \                                                     /   |
	|                                                            |
	|    \                                                 /     |
	|                                                            |
	|     \                                               /      |
	|       -- \                                     - -/        |
	|           \                                   /            |
	|                                                            |
	|             \                               /              |
	|                                                            |
	|              \                             /               |
	|                -- \                   - -/                 |
	|                    \                 /                     |
	|                      \             /                       |
	|                                                            |
	|                       \           /                        |
	|                                                            |
	|                         \- \/- //                          | 1.00
	-------------------------------------------------------------
	-20.00                                                    20.00


	-------------------------------------------------------------
	|^                                                          ^| 20.41
	| ^                                                       ^  |
	|  ^                                                     ^   |
	|                                                            |
	|    ^                                                 ^     |
	|                                                            |
	|     ^                                               ^      |
	|       ^^ ^                                     ^ ^^        |
	|           ^                                   ^            |
	|                                                            |
	|             ^                               ^              |
	|                                                            |
	|              ^                             ^               |
	|                ^^ ^                   ^ ^^                 |
	|                    ^                 ^                     |
	|                      ^             ^                       |
	|                                                            |
	|                       ^           ^                        |
	|                                                            |
	|                         ^^ ^^^ ^^                          | 1.00
	-------------------------------------------------------------
	-20.00                                                    20.00

=#

for i=-100:100 @printf("i= %9.4f  ff= %9.4f\n",i,ff_easy(i)) end 

test_ctrl = Pikaia.set_ctrl_default(123456)
tic(); r=Pikaia.pikaia(ff_rescaling, 1, test_ctrl) ; toc();r
Pikaia.result_rescaling(r, 0.,1., -20., 20.)


# TODO Iter = 1 | Mean = -10.30292 | Best = -1.106484
 @profile Pikaia.pikaia(ff_rescaling, 1, test_ctrl)
 Profile.print()
#=
311 client.jl; _start; line: 429
        311 client.jl; run_repl; line: 186
           311 client.jl; eval_user_input; line: 109
              311 profile.jl; anonymous; line: 14
                 34 ...l/src/pikaia.jl; pikaia; line: 1015
                  2  ...l/src/pikaia.jl; select; line: 478
                  2  ...l/src/pikaia.jl; select; line: 481
                  29 ...l/src/pikaia.jl; select; line: 482
                 24 ...l/src/pikaia.jl; pikaia; line: 1018
                  1  ...l/src/pikaia.jl; select; line: 478
                  4  ...l/src/pikaia.jl; select; line: 481
                  19 ...l/src/pikaia.jl; select; line: 482
                 91 ...l/src/pikaia.jl; pikaia; line: 1025
                  3  ...l/src/pikaia.jl; encode; line: 656
                  9  ...l/src/pikaia.jl; encode; line: 658
                  1  ...l/src/pikaia.jl; encode; line: 662
                  4  ...l/src/pikaia.jl; encode; line: 664
                   2 array.jl; reverse; line: 1148
                   1 range.jl; vcat; line: 343
                   1 range.jl; vcat; line: 346
                  1  ...l/src/pikaia.jl; encode; line: 665
                  2  ...l/src/pikaia.jl; encode; line: 666
                   2 ...l/src/pikaia.jl; fortran_int; line: 78
                  1  ...l/src/pikaia.jl; fortran_int; line: 78
                  1  array.jl; getindex; line: 311
                   1 abstractarray.jl; checkbounds; line: 95
                    1 abstractarray.jl; checkbounds; line: 74
                  69 array.jl; getindex; line: 312
                 49 ...l/src/pikaia.jl; pikaia; line: 1026
                  7  ...l/src/pikaia.jl; encode; line: 656
                  4  ...l/src/pikaia.jl; encode; line: 658
                  1  ...l/src/pikaia.jl; encode; line: 662
                  22 ...l/src/pikaia.jl; encode; line: 664
                   6  array.jl; reverse; line: 1148
                   15 range.jl; vcat; line: 343
                  1  ...l/src/pikaia.jl; encode; line: 665
                  6  ...l/src/pikaia.jl; encode; line: 666
                   2 ...l/src/pikaia.jl; fortran_int; line: 78
                  4  array.jl; getindex; line: 311
                   4 abstractarray.jl; checkbounds; line: 95
                    3 abstractarray.jl; checkbounds; line: 74
                    1 range.jl; maximum; line: 101
                  4  array.jl; getindex; line: 312
                 6  ...l/src/pikaia.jl; pikaia; line: 1029
                  1 ...l/src/pikaia.jl; cross!; line: 761
                   1 ...l/src/pikaia.jl; fortran_int; line: 78
                  2 ...l/src/pikaia.jl; cross!; line: 776
                   1 range.jl; vcat; line: 343
                  1 ...l/src/pikaia.jl; cross!; line: 777
                  1 ...l/src/pikaia.jl; cross!; line: 783
                  1 ...l/src/pikaia.jl; dbg; line: 56
                 1  ...l/src/pikaia.jl; pikaia; line: 1030
                  1 ...l/src/pikaia.jl; mutate!; line: 870
                 2  ...l/src/pikaia.jl; pikaia; line: 1031
                  2 ...l/src/pikaia.jl; mutate!; line: 870
                 22 ...l/src/pikaia.jl; pikaia; line: 1034
                  2  ...l/src/pikaia.jl; decode; line: 684
                   2 array.jl; getindex; line: 151
                  11 ...l/src/pikaia.jl; decode; line: 689
                  1  ...l/src/pikaia.jl; decode; line: 695
                  4  ...l/src/pikaia.jl; decode; line: 698
                   3 array.jl; push!; line: 660
                   1 array.jl; push!; line: 661
                  1  array.jl; setindex!; line: 495
                  3  array.jl; setindex!; line: 496
                   3 abstractarray.jl; checkbounds; line: 95
                    3 abstractarray.jl; checkbounds; line: 74
                 20 ...l/src/pikaia.jl; pikaia; line: 1035
                  3 ...l/src/pikaia.jl; decode; line: 684
                   3 array.jl; getindex; line: 151
                  9 ...l/src/pikaia.jl; decode; line: 689
                  2 ...l/src/pikaia.jl; decode; line: 698
                   2 array.jl; push!; line: 660
                  6 array.jl; setindex!; line: 496
                   6 abstractarray.jl; checkbounds; line: 95
                    6 abstractarray.jl; checkbounds; line: 74
                 24 ...l/src/pikaia.jl; pikaia; line: 1043
                  24 ...l/src/pikaia.jl; generational_replacement; line: 349
                   17 array.jl; copy!; line: 51
                    17 array.jl; unsafe_copy!; line: 38
                     17 array.jl; unsafe_copy!; line: 31
                 36 ...l/src/pikaia.jl; pikaia; line: 1054
                  1  ...l/src/pikaia.jl; new_pop!; line: 597
                  1  ...l/src/pikaia.jl; new_pop!; line: 612
                  32 ...l/src/pikaia.jl; new_pop!; line: 618
                   15 array.jl; getindex; line: 311
                    14 abstractarray.jl; checkbounds; line: 95
                     13 abstractarray.jl; checkbounds; line: 74
                       1  int.jl; <; line: 228
                   7  array.jl; getindex; line: 312
                     1 array.jl; setindex!; line: 412
                     1 none; ff_rescaling; line: 3
                     1 none; ff_rescaling; line: 6
                  2  ...l/src/pikaia.jl; new_pop!; line: 622
                   2 ...l/src/pikaia.jl; rank_pop; line: 526
                      2 sort.jl; sort!; line: 444
                       1 sort.jl; fpsort!; line: 426
                        1 sort.jl; nans2right!; line: 403
                       1 sort.jl; fpsort!; line: 438
                        1 sort.jl; sort!; line: 281
                 1  ...l/src/pikaia.jl; pikaia; line: 1066
                    1 ...l/src/pikaia.jl; report; line: 141
                       1 stream.jl; write; line: 679
=#
Profile.print(format=:flat)

test_ctrl = Pikaia.set_ctrl_default(sort(123456))
Pikaia.pikaia(ff_rescaling, 1, test_ctrl)

function rastrigin_rescaling(x)
#   appropriately rescaling x
    sx = x[1]
    rx = Pikaia.rescaling(sx,0.,1., -5., 5.) 
#   search min function
    return -1.*TestFunctions.rastrigin([rx])
end # rastrigin_rescaling




# *********************************************************************
function create_vector(minx, maxx, count_elements)
# =====================================================================
# creation of a vector of real numbers of a given size in the certain range
# =====================================================================
    if maxx < minx
        error("create_vector:: maxx < minx")
    end

    step = (maxx-minx)/count_elements
    xs   = Float64[]
    x  = 0.

    while x <= 1.
        push!(xs,x)
        x = x+step
    end

    xs
end
using ASCIIPlots
scatterplot(xs, map(rastrigin_rescaling, xs))

#=
    
	-------------------------------------------------------------
	|                             ^                              | -0.00
	|                                                            |
	|                        ^         ^                         |
	|                                                            |
	|                                                            |
	|                  ^                     ^                   |
	|                                                            |
	|                                                            |
	|          ^     ^                         ^     ^           |
	|     ^       ^       ^               ^       ^       ^      |
	|                          ^     ^                           |
	|                                                            |
	|^                                                          ^|
	|                                                            |
	|                                                            |
	|        ^                                         ^         |
	|                                                            |
	|                                                            |
	|                                                            |
	|  ^                                                     ^   | -40.26
	-------------------------------------------------------------
	0.00                                                    1.00

=#

test_ctrl = Pikaia.set_ctrl_default(423456)
res=Pikaia.pikaia(rastrigin_rescaling, 1, test_ctrl)
Pikaia.result_rescaling(res, 0.,1.,-5.,5.)

# (-0.012221665193063558,0.02961916178932178,0)

