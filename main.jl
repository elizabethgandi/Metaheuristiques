# ===========================================================================#
# Compilant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra
using Random


include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("heuristiqueDeConstruction.jl")
include("simpleDescent.jl")
include("HAmelioration.jl")

# ===========================================================================#

# Loading a SPP instance
println("\nLoading...")
fname = "Data/didactic.dat"
#fname = "Data/pb_2000rnd0100.dat"
#fname = "Data/pb_1000rnd0700.dat"
#fname = "Data/pb_500rnd1300.dat"
C, A = loadSPP(fname)

# ============================================================================
#PARTIE DM1===================================================================
println("\nDM1...\n")

# PARTIE CONSTRUCTION---------------------------------------------------------
println("\nSolving with construction Glouton...\n")

@time x, zBest = gloutonConstruction(C,A)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest)

# PARTIE AMELIORATION: PLUS PROFONDE DESCENTE---------------------------------
print("\nSolving with Amelioration Glouton...")

@time xAmelioration, zAmelioration = gloutonAmelioration(C, A, x, zBest)
println("x[i]=1 en i ∈ ", findall(x->x==1, xAmelioration))
println("z(x) = ", zAmelioration)

# ============================================================================
# ============================================================================

#PARTIE DM2===================================================================
println("\nDM2...\n")

#PARTIE METAHEURISTIQUES: GRASP----------------------------------------------
println("\nSolving with GRASP only...\n")

nbIterations = 10
alpha        = 0.700

@time xGRASP, zGRASP = GRASP(C, A, alpha, nbIterations)
println("x[i]=1 en i ∈ ", findall(x->x==1, xGRASP))
println("z(x) = ", zGRASP)

#=
# PARTIE METAHEURISTIQUES: GRASP AVEC DESTROY AND REPAIR----------------------
println("\nSolving with GRASP with destroy and repair...\n")

nbIterationsDR = 10
alphaDR        = 0.700

=#

#=
println("\nSolving with Destroy and repear...\n")

x, zBest = destroyAndRepear(Ctemp, A)
println("x[i]=1 en i ∈ ", findall(x->x==1, x))
println("z(x) = ", zBest3)


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
