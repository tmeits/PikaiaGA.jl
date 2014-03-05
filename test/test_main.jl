#

using Base.Test
using Pikaia

a = rand(1 : 50)

@test sort(a) == a[Pikaia.rqsort(length(a), a, [1:50])]

@test Pikaia.setctl([100, 500, 7, .85, 0, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 5
@test Pikaia.setctl([100, 500, 7, .85, 7, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 5
@test Pikaia.setctl([100, 500, 7, .85, 2, .005, .0005, .25, 1, 1, 1, 0], 2)[1] == 0


