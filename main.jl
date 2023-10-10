# =========================================================================== #
# Compilant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("heuristiqueDeConstruction.jl")

# =========================================================================== #

# Loading a SPP instance
println("\nLoading...")
#fname = "Data/didactic.dat"
#
fname = "Data/pb_500rnd0100.dat"
C, A = loadSPP(fname)

println("\nSolving with Glouton...")
Ctemp = copy(C)
@time x, z = glouton(Ctemp,A)

@show z
@show x

# Solving a SPP instance with GLPK
#=println("\nSolving with GLPK...")
solverSelected = GLPK.Optimizer
spp = setSPP(C, A)

set_optimizer(spp, solverSelected)
optimize!(spp)

# Displaying the results
println("z = ", objective_value(spp))
print("x = "); println(value.(spp[:x]))

=#


# =========================================================================== #

# Collecting the names of instances to solve
#=println("\nCollecting...")
target = "Data"
fnames = getfname(target)=#


println("\nThat's all folks !")
