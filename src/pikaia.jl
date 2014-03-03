
module Pikaia

importall Base

export pikaia, rqsort

function pikaia(ff:: Function, n:: Int, ctrl:: Array{Float64, 1})
# Optimization (maximization) of user-supplied "fitness" function
# ff  over n-dimensional parameter space  x  using a basic genetic
# algorithm method.

# Version 0.0.1   [ 2014 February 21 ]

# Init output:
x = Array{Float64, n}
f :: Float64 = 0.
status :: Int = 0

x = rand(n)
println(n, x, ff(n, x))
# Return in typle
return (x, f, status)
end # pikaia

function rqsort(n:: Int, a:: Array{Float64, 1}, p:: Array{Int, 1})
# Return integer array p which indexes array a in increasing order.
# Array a is not disturbed.
    p=sortperm(a)
    return p
end

end # Pikaia
#
# #3F3434 #DFE8E2
# http://habrahabr.ru/post/161009/
# https://github.com/johnmyleswhite/HopfieldNets.jl

