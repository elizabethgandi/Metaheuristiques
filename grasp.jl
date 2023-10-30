function GRASP(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, alpha::Float64, nbIterations::Int64)

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation ----------------------------------------------------------
    n::Int64               = size(A,2)       # n nombre de variables
    m::Int64               = size(A,1)       # m nombre de contraintes

    index::Vector{Int64}   = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}     = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64               = 0               # z la valeur de la fonction objective

    randomCandidate::Int64 = 0
    zMeilleur              = 0
    xMeilleur              = zeros(Int64,n)


    # Partie GRASP construction -----------------------------------------------

    # 1) REDUCTION DE L'INSTANCE SUR VARIABLES NON CONTRAINTES ----------------

    # Elimine de l'instance toutes les variables concernees par aucune contrainte
    variablesConcernees = findall(isequal(0), vec(sum(A, dims=1)))
    for j in variablesConcernees
        sol[j] = 1      # maj solution (partielle) en fixant a 1 la variable non contrainte
        z = z + C[j]    # maj valeur de la solution (partielle)
    end

    # Supprime les colonnes correspondant aux variables fixees
    index = index[setdiff(1:end, variablesConcernees)]
    C = C[setdiff(1:end, variablesConcernees)]
    A = A[:, setdiff(1:end, variablesConcernees)]

    # 2) CALCUL DES UTILITES --------------------------------------------------

    candidates::Vector{Float64} = C ./ vec(sum(A, dims = 1))

    if verbose
        println("ivar   : ", index)
        #println("C      : ", C)
        #println("A      : ", A)
        println("U      : ", candidates)
        println(" ")
    end

    Ccopy          = copy(C)
    Acopy          = copy(A)
    candidatesCopy = copy(candidates)
    indexCopy      = copy(index)
    zCopy          = z
    solCopy        = copy(sol)

    for i in 1:nbIterations
        print("\nItération n° $i: ")

        C          = copy(Ccopy) 
        A          = copy(Acopy) 
        candidates = copy(candidatesCopy)
        index      = copy(indexCopy)
        z          = zCopy
        sol        = copy(solCopy)

        while (size(A,1) != 0) && (size(C,1) != 0) 

            print("x")

            limite = candidates[argmin(candidates)] + (alpha)*(candidates[argmax(candidates)]-candidates[argmin(candidates)])
            RCL    = findall(x-> (x >= limite), candidates)

            # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE -----------

            # Selection au hazard de l'indice d'un des candidats 
            randomCandidate = RCL[rand(1:size(RCL, 1))]

            # Mise à jour de la solution avec le candidat selectionne
            sol[index[randomCandidate]] = 1
            # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
            z = z + C[randomCandidate]

            # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE ----------

            # Identification des contraintes a supprimer du fait de la variable selectionnee
            lignestemp = findall(isequal(1), A[:,randomCandidate])

            # identifie toutes les colonnes qui doivent etre supprimées suite au candidat selectionne
            colonnetemp=(Int64)[]
            for i in lignestemp
                # scrute contrainte par contrainte les coefficients de A de valeur 1 (et elimine les eventuels doublons)
                colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
            end

            # On supprime dans les structures index, A et C les valeurs correspondantes aux colonnes supprimées
            index      = index[setdiff(1:end, colonnetemp)]      # reduction du vecteur des indices des variables
            A          = A[:, setdiff(1:end, colonnetemp)]       # reduction des colonnes de la matrice des contraines
            C          = C[setdiff(1:end, colonnetemp)]          # reduction du vecteur de coefficients de la fonction objectif
            candidates = candidates[setdiff(1:end, colonnetemp)] # reduction du vecteur des utilites

            # On supprime dans la matrice A les lignes correspondantes aux contraintes impliquees par la variable selectionnee
            A          = A[setdiff(1:end, lignestemp), :]        # reduction des lignes de la matrice des contraines
        
            if verbose
                println("jselec : ", randomCandidate)
                println("ivar   : ", index)
                #println("C      : ", C)
                #println("A      : ", A)
                println("U      : ", candidates)
                println("-------- ")
            end
            nbIterations -= 1

        end

        if(verbose)
            print("zGRASP = $z \n")
        end

        # Partie GRASP amélioration -------------------------------------------

        xA, zA = gloutonAmelioration(C_entree, A_entree, sol, z)

        if (zA > zMeilleur)
            zMeilleur = zA
            xMeilleur = xA
        end

        if(verbose)
            print("zAmelioration = $zA \n")
        end

    end

    if(verbose)
        println(" ")
    end
    return sol, zMeilleur
end


function grasp_DR(C, A, alpha, nbIterationGrasp, nbIterationDR)

    zconstruction = zeros(Int64,nbIterationGrasp)
    zamelioration = zeros(Int64,nbIterationGrasp)
    zbest         = zeros(Int64,nbIterationGrasp)
    zbetter       = 0
    zMax          = 0

    zDeconstruit  = zeros(Int64,nbIterationDR)
    zReconstruit  = zeros(Int64,nbIterationDR)
    xReconstruit::Vector{Int64} = zeros(Int64, size(C,1))

    # Iterations GRASP --------------------------------------------------------
    for i=1:nbIterationGrasp

        @printf("\niteration GRASP n° %s  \n", i)

        # Construction grasp --------------------------------------------------
        print("  construction grasp : ")
        x, zconstruction[i] = graspSPP(C, A, alpha, nbIterationGrasp)
        println(" zGraspC = ", zconstruction[i])

        # Amelioration grasp --------------------------------------------------
        print("  amelioration grasp : ")
        x, zamelioration[i] = gloutonAmelioration(C, A, x, zconstruction[i]) 
        println(" zGraspA = ", zamelioration[i])

        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
        
        # Iterations Deconstruction/Reconstruction ----------------------------
        for j=1:nbIterationDR

            @printf("    iteration DR n° %s  \n", j)

            # Deconstruction --------------------------------------------------
            print("    deconstruction   : ")
            xDeconstruit, zDeconstruit[j] = destruction(C, x, zamelioration[i])
            println(" zDecons = ", zDeconstruit[j])

            # Reconstruction --------------------------------------------------
            println("    reconstruction   :   ")
            xReconstruit, zReconstruit[j] = reconstructionDR(C, A, xDeconstruit, zDeconstruit[j])
            #xReconstruction, zReconstruction = GreedyAmeliorationVND(C, A, xDeconstruit, zDeconstruit)
            println("                        zRecons = ", zReconstruit[j])

            # Mise a jour si besoin de la meilleure solution rencontree -------
            zbetter = max(zbetter, zReconstruit[j])
            zbest[i] = zbetter
        end
        
    end

    println("\n zbetter : ", zbetter)

    return xReconstruit, zbetter
end


# --------------------------------------------------------------------------- #

function selectionRCL(u, alpha)
    verbose = false

    epsilon = 10^-10
    umin = minimum(u)
    umax = maximum(u)
    ulimit = umin + alpha * (umax - umin)
    rcl = (Int64)[]
    for j=1:length(u)
        if u[j] >= ulimit - epsilon
            push!(rcl,j)
        end
    end
    jselect = rcl[rand(1:length(rcl))]

    if verbose
        @show u
        @printf(" %5.2f  %5.2f  %5.2f \n", umin, ulimit, umax)
        @show rcl
        @printf(" %3d %5.2f \n", jselect, u[jselect])
    end
    return jselect
end

# --------------------------------------------------------------------------- #

function graspSPPconstruction(C_in, A_in, alpha)
    C = deepcopy(C_in)
    A = deepcopy(A_in)
    verbose = false

    z::Int64 = 0  # valeur de la solution
    x = zeros(Int64, length(C)) # variables de la solution a priori toutes a zero
    posVar = collect(1:1:length(C)) # vecteur d'indices des variables
    ligneSaturee = Int64[]

    #println("Construction gloutonne evaluee d'une solution admissible")

    # Elimine de l'instance toutes les variables concernees par aucune contrainte
    variablesConcernees = findall(isequal(0), vec(sum(A, dims=1)))
    for j in variablesConcernees
        x[j] = 1 # maj solution (partielle) en fixant a 1 la variable non contrainte
        z = z + C[j] # maj valeur de la solution (partielle)
    end
    # supprime les colonnes correspondant aux variables fixees
    posVar = posVar[setdiff(1:end, variablesConcernees)]
    C = C[setdiff(1:end, variablesConcernees)]
    A = A[:, setdiff(1:end, variablesConcernees)]

    while (size(A,1) != 0) && (size(C,1) != 0)
        utilite = C ./ vec(sum(A, dims=1))

# !!!!!!!!!!!!! DEBUT : Passage qui diffère avec le tout glouton !!!!!!!!!!!!!
        j_selec = selectionRCL( utilite , alpha )
# !!!!!!!!!!!!! FIN : Passage qui diffère avec le tout glouton !!!!!!!!!!!!!

        x[posVar[j_selec]] = 1 # maj solution (partielle)
        z = z + C[j_selec] # maj valeur de la solution (partielle)

        if verbose
            println("ivar   : ", posVar)
            println("C      : ", C)
            println("A      : ", A)
            println("U      : ", utilite)
            println("jselec : ",j_selec)
            println("-------- ")
        end

       # Elimine toutes les variables fixees et reduit l'instance
       ligneSaturee = findall(isequal(1), A[:,j_selec]) # ligne de A baree
       variablesConcernees = Int64[]
        for i in ligneSaturee # parcours de la colonne selectionnee
            variablesConcernees = union(variablesConcernees,findall(isequal(1), A[i,:])) # lignes de A concernees
        end
        # supprime les colonnes correspondant aux variables fixees
        variablesConcernees = unique(variablesConcernees) # elimine les indices doublons
        posVar = posVar[setdiff(1:end, variablesConcernees)]
        C = C[setdiff(1:end, variablesConcernees)]
        A = A[:, setdiff(1:end, variablesConcernees)]

        # Elimine toutes les contraintes saturees et reduit l'instance
        A = A[setdiff(1:end, ligneSaturee), :]

        # Elimine toutes les contraintes de A contenant que des zeros
        contraintesConcernees = Int64[]
        for ligne in (1:size(A,1))
            if findfirst(isequal(1), A[ligne,:]) == nothing
                contraintesConcernees = union(contraintesConcernees,ligne)
            end
        end
        A = A[setdiff(1:end, contraintesConcernees), :]
    end
    return x, z
end


# --------------------------------------------------------------------------- #

function graspSPP(C, A, alpha, nbIterationGrasp)

    zconstruction = zeros(Int64,nbIterationGrasp)
    zamelioration = zeros(Int64,nbIterationGrasp)
    zbest         = zeros(Int64,nbIterationGrasp)
    zbetter       = 0

    for i=1:nbIterationGrasp
        x, zconstruction[i] = graspSPPconstruction(C, A, alpha)
        x, zamelioration[i] = GreedyAmelioration(C, A, x, zconstruction[i]) # rand(0:10) + zconstruction[i] # livrable du DM2
        @printf("z(xGrasp) = %d | iteration = %s \n",zamelioration[i], i)
        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
    end
    return zconstruction, zamelioration, zbest
end

# --------------------------------------------------------------------------- #

function graspSPP_DR(C, A, alpha, nbIterationGrasp, nbIterationDR)

    zconstruction = zeros(Int64,nbIterationGrasp)
    zamelioration = zeros(Int64,nbIterationGrasp)
    zbest         = zeros(Int64,nbIterationGrasp)
    zbetter       = 0
    zMax          = 0

    zDeconstruit  = zeros(Int64,nbIterationDR)
    zReconstruit  = zeros(Int64,nbIterationDR)
    xReconstruit::Vector{Int64} = zeros(Int64, size(C,1))

    # Iterations GRASP --------------------------------------------------------
    for i=1:nbIterationGrasp

        @printf("\niteration GRASP n° %s  \n", i)

        # Construction grasp --------------------------------------------------
        print("  construction grasp : ")
        x, zconstruction[i] = graspSPPconstruction(C, A, alpha)
        #x, zconstruction[i] = GRASP(C, A, alpha, nbIterationGrasp)
        println(" zGraspC = ", zconstruction[i])

        # Amelioration grasp --------------------------------------------------
        print("  amelioration grasp : ")
        #x, zamelioration[i] = GreedyAmelioration(C, A, x, zconstruction[i]) 
        x, zamelioration[i] = gloutonAmelioration(C, A, x, zconstruction[i]) 
        println(" zGraspA = ", zamelioration[i])

        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
        
        # Iterations Deconstruction/Reconstruction ----------------------------
        for j=1:nbIterationDR

            @printf("    iteration DR n° %s  \n", j)

            # Deconstruction --------------------------------------------------
            print("    deconstruction   : ")
            xDeconstruit, zDeconstruit[j] = destruction(C, x, zamelioration[i])
            println(" zDecons = ", zDeconstruit[j])

            # Reconstruction --------------------------------------------------
            println("    reconstruction   :   ")
            #xReconstruit, zReconstruit[j] = GreedyReconstruction(C, A, xDeconstruit, zDeconstruit[j])
            xReconstruit, zReconstruit[j] = reconstructionDR(C, A, xDeconstruit, zDeconstruit[j])
            
            #xReconstruction, zReconstruction = GreedyAmeliorationVND(C, A, xDeconstruit, zDeconstruit)
            println("                        zRecons = ", zReconstruit[j])

            # Mise a jour si besoin de la meilleure solution rencontree -------
            zbetter = max(zbetter, zReconstruit[j])
            zbest[i] = zbetter
        end
        
    end

    println("\n zbetter : ", zbetter)

    return xReconstruit, zbetter
end