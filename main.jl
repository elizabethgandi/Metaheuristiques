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
fname = "Data/pb_1000rnd0100.dat"
C, A = loadSPP(fname)


println("\nSolving with Glouton...\n")
Ctemp = copy(C)

@time x, zBest = glouton(Ctemp,A)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest)

println("\nSolving with GRASP...\n")

@time x, zBest = GRASP(Ctemp, A)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest)


#=println("\nUsing simple descent to upgrade the solution... may take time")
xnew = simpleDescent(copy(C), copy(A), copy(x))
@show xnew

#println("\nSolving with Glouton...")
Ctemp = copy(C)

#@time x, zOpt = glouton(Ctemp,A)
x, zOpt = glouton(Ctemp,A)

#@show zOpt
#@show x


#=Solving a SPP instance with GLPK
println("\nSolving with GLPK...")
solverSelected = GLPK.Optimizer
spp = setSPP(C, A)

set_optimizer(spp, solverSelected)
@time optimize!(spp)

# Displaying the results
println("z = ", objective_value(spp))
print("x = "); println(value.(spp[:x]))
=#


#println("\nUsing simple descent to upgrade the solution... may take time")
#@time xnew = simpleDescent(copy(C), copy(A), copy(x))
#@show xnew


println("\nUsing another simple descent to upgrade the solution")
@time xnew2 = simpleDescent2(copy(C), copy(A), copy(x))
@show xnew2
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

# =========================================================================== 

#Collecting the names of instances to solve

println("\nCollecting...")
target = "Data"
fnames = getfname(target)
=#

println("\nThat's all folks !")
