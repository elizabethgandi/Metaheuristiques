# =========================================================================== #
# Compilant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("heuristiqueDeConstruction.jl")
include("simpleDescent.jl")

# =========================================================================== #

# Loading a SPP instance
println("\nLoading...")
fname = "Data/didactic.dat"
C, A = loadSPP(fname)
@show C
@show A

println("\nSolving with Glouton...")
Ctemp = copy(C)
x = glouton(Ctemp,A)
@show x

println("\nUsing simple descent to upgrade the solution")
xnew = simpleDescent(copy(C), copy(A), copy(x))
@show xnew

# Solving a SPP instance with GLPK
println("\nSolving with GLPK...")
solverSelected = GLPK.Optimizer
spp = setSPP(C, A)

set_optimizer(spp, solverSelected)
optimize!(spp)

# Displaying the results
println("z = ", objective_value(spp))
print("x = "); println(value.(spp[:x]))


# =========================================================================== #

# Collecting the names of instances to solve
#=println("\nCollecting...")
target = "Data"
fnames = getfname(target)=#




println("\nThat's all folks !")
