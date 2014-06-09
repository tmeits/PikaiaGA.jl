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

function f(x)
    return abs(x) + cos(x)
end    


