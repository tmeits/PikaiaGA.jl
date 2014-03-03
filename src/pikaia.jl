
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

function setctl(ctrl:: Array{Float64, 1})
# Set control variables and flags from input and defaults
const DEFAULT = [100, 500, 5, .85, 2, .005, .0005, .25, 1, 1, 1, 0]
    for i=1:length(ctrl)
        if ctrl[i] < 0.
            ctrl[i] = DEFAULT[i]
        end
    end
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
    setctl(ctrl)

#println(n, x, ff(n, x))

# Return in typle
return (x, f, status)
end # pikaia

end # Pikaia
#
# #3F3434 #DFE8E2
# http://habrahabr.ru/post/161009/
# https://github.com/johnmyleswhite/HopfieldNets.jl

