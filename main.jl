# =========================================================================== #
# Compilant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra
#using (SparseArrays)

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("heuristiqueDeConstruction.jl")
include("simpleDescent.jl")

# =========================================================================== #

# Loading a SPP instance
println("\nLoading...")
#fname = "Data/didactic.dat"
#fname = "Data/pb_100rnd0100.dat"
#fname = "Data/pb_100rnd0200.dat"
fname = "Data/pb_500rnd0500.dat"
#fname = "Data/pb_2000rnd0100.dat"


C, A = loadSPP(fname)

#=
C = [1, 1, 1, 5]
A = [ 1 0 0 1;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1]
=#

println("\nSolving with Glouton...\n")
Ctemp = copy(C)

@time x, zBest = glouton(Ctemp,A)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest)

println("\nImproving the solution with a simple descent...\n")
@time x, zBest = simpleDescent(C,A,x,zBest)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest)


#=

#Solving a SPP instance with GLPK
println("\nSolving with GLPK...")
solverSelected = GLPK.Optimizer
spp = setSPP(C, A)

set_optimizer(spp, solverSelected)
@time optimize!(spp)

# Displaying the results
println("z = ", objective_value(spp))
print("x = "); println(value.(spp[:x]))

=#

#= Test fun 
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

=#

# =========================================================================== 

#=

#Collecting the names of instances to solve

println("\nCollecting...")
target = "Data"
fnames = getfname(target)


println("\nThat's all folks !")

=#