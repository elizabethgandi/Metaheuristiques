
function grasp_DR(C, A, alpha, nbIterationGrasp, nbIterationDR)

    zconstruction = 0 
    zamelioration = 0 
    zbetter       = 0
    zMax          = 0

    zDeconstruit  = 0 
    zReconstruit  = 0 
    xReconstruit::Vector{Int64} = zeros(Int64, size(C,1))

    # Iterations GRASP --------------------------------------------------------
    for i=1:nbIterationGrasp

        @printf("\niteration GRASP n° %s  \n", i)

        # Construction grasp --------------------------------------------------
        print("  construction grasp : ")
        xconstruction, zconstruction = GRASPconstruction(C, A, alpha)
        println(" zGraspC = ", zconstruction)

        # Amelioration grasp --------------------------------------------------
        print("  amelioration grasp : ")
        xamelioration, zamelioration = gloutonAmelioration(C, A, xconstruction, zconstruction) 
        println(" zGraspA = ", zamelioration)

        zbetter = max(zbetter, zamelioration)
        
        # Iterations Deconstruction/Reconstruction ----------------------------
        for j=1:nbIterationDR

            @printf("    iteration DR n° %s  \n", j)

            # Deconstruction --------------------------------------------------
            print("    deconstruction   : ")
            xDeconstruit, zDeconstruit = destruction(C, xamelioration, zamelioration)
            println(" zDecons = ", zDeconstruit)

            # Reconstruction --------------------------------------------------
            println("    reconstruction   :   ")
            xReconstruit, zReconstruit = reconstructionDR(C, A, xDeconstruit, zDeconstruit)
            println("                        zRecons = ", zReconstruit)

            # Mise a jour si besoin de la meilleure solution rencontree -------
            zbetter = max(zbetter, zReconstruit)
        end
    end

    println("\n zbetter : ", zbetter)

    return xReconstruit, zbetter
end



function GRASP(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, alpha, nbIterations)

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation
    n::Int64             = size(A,2)       # n nombre de variables
    m::Int64             = size(A,1)       # m nombre de contraintes

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    randomCandidate::Int64 = 0
    zMeilleur              = 0
    xMeilleur              = zeros(Int64,n)


    # 1) REDUCTION DE L'INSTANCE SUR VARIABLES NON CONTRAINTES ----------------

    # Elimine de l'instance toutes les variables concernees par aucune contrainte
    variablesConcernees = findall(isequal(0), vec(sum(A, dims=1)))
    for j in variablesConcernees
        sol[j] = 1 # maj solution (partielle) en fixant a 1 la variable non contrainte
        z = z + C[j] # maj valeur de la solution (partielle)
    end

    # supprime les colonnes correspondant aux variables fixees
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

            limite = candidates[argmin(candidates)] + (alpha)*(candidates[argmax(candidates)]-candidates[argmin(candidates)])
            RCL = findall(x-> (x >= limite), candidates)

            # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

            # Selection au hazard de l'indice d'un des candidats 
            randomCandidate = RCL[rand(1:size(RCL, 1))]

            # Mise à jour de la solution avec le candidat selectionne
            sol[index[randomCandidate]] = 1
            # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
            z = z + C[randomCandidate]

            # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

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
            #nbIterations -= 1

        end

        print("zGRASP = $z \n")

        xA, zA = gloutonAmelioration(C_entree, A_entree, sol, z)

        if (zA > zMeilleur)
            zMeilleur = zA
            xMeilleur = xA
        end

        print("zAmelioration = $zA \n")

    end

    println(" ")
    return sol, zMeilleur
end



function GRASPconstruction(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, alpha)

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation
    n::Int64             = size(A,2)       # n nombre de variables
    m::Int64             = size(A,1)       # m nombre de contraintes

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    randomCandidate::Int64 = 0


    # 1) REDUCTION DE L'INSTANCE SUR VARIABLES NON CONTRAINTES ----------------

    # Elimine de l'instance toutes les variables concernees par aucune contrainte
    variablesConcernees = findall(isequal(0), vec(sum(A, dims=1)))
    for j in variablesConcernees
        sol[j] = 1 # maj solution (partielle) en fixant a 1 la variable non contrainte
        z = z + C[j] # maj valeur de la solution (partielle)
    end

    # supprime les colonnes correspondant aux variables fixees
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

    while (size(A,1) != 0) && (size(C,1) != 0) 

        limite = candidates[argmin(candidates)] + (alpha)*(candidates[argmax(candidates)]-candidates[argmin(candidates)])
        RCL = findall(x-> (x >= limite), candidates)

        # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection au hazard de l'indice d'un des candidats 
        randomCandidate = RCL[rand(1:size(RCL, 1))]

        # Mise à jour de la solution avec le candidat selectionne
        sol[index[randomCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[randomCandidate]

        # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

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

    end

    print("zGRASP construction = $z \n")

    println(" ")
    return sol, z #sol, zMeilleur
end