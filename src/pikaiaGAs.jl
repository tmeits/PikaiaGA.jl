
module pikaiaGAs

importall Base

export pikaia, seed

#  Common block to make iseed visible to rninit (and to save
type seed
    iseed::Int
end

function pikaia(ff, n, ctrl, xb, fb, status)
# Optimization (maximization) of user-supplied "fitness" function
# ff  over n-dimensional parameter space  x  using a basic genetic
# algorithm method.

# Version 0.0.1   [ 2014 February 21 ]

end

end
#
