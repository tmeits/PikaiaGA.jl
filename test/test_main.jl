#

using Base.Test
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


