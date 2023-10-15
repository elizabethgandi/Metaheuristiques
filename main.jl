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
#fname = "Data/didactic.dat"
fname = "Data/pb_2000rnd0100.dat"
C, A = loadSPP(fname)

println("\nSolving with Glouton...")
Ctemp = copy(C)

@time x, zOpt = glouton(Ctemp,A)

@show zOpt
@show x
#x = glouton(Ctemp,A)


#=println("\nUsing simple descent to upgrade the solution... may take time")
xnew = simpleDescent(copy(C), copy(A), copy(x))
@show xnew

println("\nUsing another simple descent to upgrade the solution")
xnew2 = simpleDescent2(copy(C), copy(A), copy(x))
@show xnew2

=#
#=
# Test fun 
xtest = simpleDescent(copy(C), copy(A), [0,0,1,1,0,1,0,0])
@show xtest
=#

#=
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
println("\nCollecting...")
target = "Data"
fnames = getfname(target)=#


println("\nThat's all folks !")
