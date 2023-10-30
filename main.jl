# ===========================================================================#
# Compilant julia 1.9

println(" ----- Resolution du SPP ----- ")

# Using the following packages
#using JuMP, GLPK
#using LinearAlgebra
using Random
using Printf

include("loadSPP.jl")
#include("setSPP.jl")
include("getfname.jl")
include("construction.jl")
include("destruction.jl")
include("simpleDescent.jl")
include("amelioration.jl")
include("grasp.jl")
include("reconstruction.jl")

#Random.seed!(100)

function resolution(fnames)

    resultats = []

    for instance = 1:length(fnames)

        C, A = loadSPP(string(target,"/",fnames[instance]))
        println("Instance : ", fnames[instance])
    
        # DM1 =====================================================================

        println("\nDM1 ----------------------------------------------------------------")
        println("\nConstruction gloutonne d'une solution admissible")
        start = time()
        xConstruction, zConstruction = gloutonConstruction(C, A)
        tConstruction = time()-start
        @printf("z(xInit) = %d | t = %f sec\n", zConstruction, tConstruction)

        println("\nAmelioration par recherche locale de type plus profonde descente")
        start = time()
        xAmelioration, zAmelioration = gloutonAmelioration(C, A, xConstruction, zConstruction)
        tAmelioration = time()-start
        @printf("z(xBest) = %d | t = %f sec \n",zAmelioration, tAmelioration)


        # DM2 =====================================================================

        println("\nDM2 ----------------------------------------------------------------")
        nbIterationGrasp = 10       # nombre d'iteration GRASP
        nbIterationDR    = 1        # Destroy/repair
        alpha            = 0.700    # alpha

        start = time()
        xbest, zbest = GRASP(C, A, alpha, nbIterationGrasp)
        tgraspSPP = time()-start
        @printf("\nzBestGrasp   = %d ", zbest)
        println("t GRASP    : ", trunc(tgraspSPP, digits=3), "sec")


        start = time()
        xfinal, zfinal = grasp_DR(C, A, alpha, nbIterationGrasp, nbIterationDR)
        tgraspSPP_DR = time()-start
        @printf("\nzBestGraspDR = %d ", zfinal)
        println("t GRASP_DR : ", trunc(tgraspSPP_DR, digits=3), "sec")


        # Sauvegarde les resultats pour cette instance ============================
        push!(resultats, (fnames[instance], zAmelioration, trunc(tAmelioration, digits=3), zbest, trunc(tgraspSPP, digits=3), zfinal, trunc(tgraspSPP_DR, digits=3)) )

    end

    return resultats
end


# =========================================================================== #
# =========================================================================== #
# Entry point

target = "Data"            # chemin vers le repertoire des instances

# experimente une instance :
#fnames = ["didactic.dat"]
fnames = ["pb_100rnd0100.dat"]
#fnames = ["pb_200rnd0900.dat"]
#fnames = ["pb_2000rnd0100.dat"]
#fnames = ["pb_500rnd0100.dat"]
#fnames = ["pb_500rnd0300.dat"]
#fnames = ["pb_2000rnd0800.dat"]

# experimente toutes les instances :
#fnames = getfname(target)

# collecte les resultats DM1 + DM2/GRASP + DM2/GRASP_DR
resultats = resolution(fnames)

println("\nEdition des resultats ----------------------------------------------")
nCA    = 0
nGRASP = 0
nDR    = 0
for r in 1:length(resultats)
    print(resultats[r])
    meilleur = max(resultats[r][2],resultats[r][4],resultats[r][6])
    print("      ")
    if meilleur == resultats[r][2]
        print(" C+A ")
        global nCA+=1
    end
    if meilleur == resultats[r][4]
        print(" GRASP ")
        global nGRASP+=1
    end
    if meilleur == resultats[r][6]
        print(" DR ")
        global nDR+=1
    end
    println(" ")
end

println("Nb instances : ", length(resultats), " | C+A : ", nCA, " | GRASP : ", nGRASP, " | DR : ", nDR)

println("that's all folk")
