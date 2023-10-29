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
include("latexTable.jl")

# =========================================================================== #

function main()
      #Collecting the names of instances to solve

      println("\nCollecting...")
      target = "Data"
      fnames = getfname(target)


      gl_times    = Vector{Float64}(undef, length(fnames))
      gl_results  = Vector{Int64}(undef, length(fnames))
      ds_times    = Vector{Float64}(undef, length(fnames))
      ds_results  = Vector{Int64}(undef, length(fnames))
      cpt         = 1


      for fname in fnames

            # Loading a SPP instance
            println("\nLoading " * fname * "..." )

            C, A = loadSPP("Data/" * fname)

            # Solving with glouton algorithm
            println("\nSolving with Glouton...\n")
            
            Ctemp = copy(C)

            gl_times[cpt] = @elapsed @time x, zBest = glouton(Ctemp,A)
            gl_results[cpt] = zBest

            println("x[i]=1 en i ∈ ", findall(x->x==1, x))
            println("z(x) = ", zBest)

            println("\nImproving the solution with a simple descent...\n")

            ds_times[cpt] = @elapsed @time x, zBest = simpleDescent(C,A,x,zBest)
            ds_results[cpt] = zBest 

            println("x[i]=1 en i ∈ ", findall(x->x==1, x))
            println("z(x) = ", zBest)

            cpt = cpt + 1
            
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

      end 

      println("\n Bilan :")

      dm1_latexTable(fnames, gl_times, gl_results, ds_times, ds_results) |> print


      # =========================================================================== 

      println("\nThat's all folks !")
end