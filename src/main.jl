# Sample driver program for pikaiaGAs.jl

import pikaiaGAs

# Initialize some empty vectors
ctrl = Float64[]

# First, initialize the random-number generator
srand(123455)

# Set control variables (use defaults)
for i= [1 : 12]
        # Push default values into an array
	ctrl = push!(ctrl, -1.)
end
ctrl[2]=50.

# Now call pikaia
res=pikaiaGAs.pikaia(twod, 2, ctrl)

# Print the results
println(" status: ", res[3])

function twod(n:: Int, x:: Array{Float64, 1})
# Compute sample fitness function (2-d landscape)
	const pi= 3.1415926536
	const sigma2= 0.15
	const nn= 9

	if x[1] > 1. || x[2] > 1.
		println("stop")
	end
	rr=sqrt( (0.5 - x[1])^2+ (0.5-x[2])^2) 
	return cos(rr*nn*pi)^2 *exp(-rr^2/sigma2) 
end
