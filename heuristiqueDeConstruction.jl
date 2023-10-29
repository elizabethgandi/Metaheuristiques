function gloutonConstruction(C_entree::Vector{Int64}, A_entree::Matrix{Int64})

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation
    n::Int64             = size(A,2)       # n nombre de variables
    m::Int64             = size(A,1)       # m nombre de contraintes

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    bestCandidate::Int64 = 0


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

        # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection de l'indice du meilleur candidat au regard de son utilite
        bestCandidate = argmax(candidates)
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]


        # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

        # Identification des contraintes a supprimer du fait de la variable selectionnee
        lignestemp = findall(isequal(1), A[:,bestCandidate])

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
            println("jselec : ", bestCandidate)
            println("ivar   : ", index)
            #println("C      : ", C)
            #println("A      : ", A)
            println("U      : ", candidates)
            println("-------- ")
        end
     
    end
    return sol, z
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
            nbIterations -= 1

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



function destroyAndRepear(C, A) # appel apres l'amelioration

    #INITIALISATION -------------------------------------------

    Ccopy                        = C
    Acopy                        = A
    xOptimal, zOptimal           = glouton(Ccopy, Acopy)
    xtemp::Vector{Int64}         = zeros(Int64, size(A, 2))
    ztemp::Int64                 = 0
    estVisite::Vector{Bool}      = fill(false, size(C, 1))

    while true

        # DESTRUCTION ------------------------------------------
        xD, zD, estVisite = detruireSolution(Ccopy, Acopy, estVisite)
        println("Destruction en cours...")
        #@show xD, zD, estVisite

        # RECONSTRUCTION ---------------------------------------
        #On veut fixer envoyer a construction C modifié pour qu'il nous donne une solution admissible
        xtemp, ztemp = reconstruction(xD, zD, A)

        if( ztemp > zOptimal)
            zOptimal = ztemp
            xOptimal = xtemp
        end

        #condition d'arret ---------------------------------------
        (all(estVisite)==false)&&break  
    end
    return xOptimal, zOptimal
end

function detruireSolution(C, A, estVisite)
    
    #xPrime, zPrime                   = simpleDescent(C, A, x) #fonctionne pas
    
    #xPrime = [0, 0, 0, 1, 0, 1, 1, 0, 0]
    #zPrime = 30

    xPrime, zPrime                   = glouton(C, A)

    L                                = findall(isequal(1), xPrime)
    shuffle!(L)

    for i in eachindex(L)
        estVisite[i] = true
    end
    
    variablesChoisies::Int64         = rand(1:size(L, 1))
    nbDeVarChoisies::Int64           = length(variablesChoisies)
    nbElemDetruits                   = rand(1:nbDeVarChoisies)
    ListeDetruits                    = Vector(1:nbElemDetruits)

    for i in eachindex(ListeDetruits)
        C[i]         = 0
        zPrime       = zPrime - C[i] #z =z-x5 
        estVisite[i] = true
    end

    return C, zPrime, estVisite

end


function reconstruction(C::Vector{Int64}, zD::Int64, A::Matrix{Int64})

    println("Reconstruction en cours...")

    verbose::Bool = false

    # Initialisation
    n::Int64                = size(A,2)       # n nombre de variables
    index::Vector{Int64}    = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}      = zeros(Int64,n)  # Vecteur de base de la solution
    lignesASupprimerDAvance = 0
    z::Int64                = 0               # z la valeur de la fonction objective
    bestCandidate::Int64    = 0

    # 0) CALCUL DES UTILITES --------------------------------------------------

    candidates::Vector{Float64} = C ./ vec(sum(A, dims = 1))

    if verbose
        println("ivar   : ", index)
        #println("C      : ", C)
        #println("A      : ", A)
        println("U      : ", candidates)
        println("a      : ", alpha)
        println(" ")
    end
    

    # PRE-TRAITEMENT --------------------------------------------------------
    for j in eachindex(C)
        lignesASupprimerDAvance = findall(isequal(1), A[:,j])
    end

    colonnetemp=(Int64)[]
    for i in lignesASupprimerDAvance
         colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
    end
    index      = index[setdiff(1:end, colonnetemp)]      
    A          = A[:, setdiff(1:end, colonnetemp)]       
    C          = C[setdiff(1:end, colonnetemp)]          
    candidates = candidates[setdiff(1:end, colonnetemp)] 
    A          = A[setdiff(1:end, lignesASupprimerDAvance), :]        
    #

    while (size(A,1) != 0) && (size(C,1) != 0) 

        # 1) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection de l'indice du meilleur candidat au regard de son utilite
        bestCandidate = argmax(candidates)

        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        zD = z + C[bestCandidate]

        # 2) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

        # Identification des contraintes a supprimer du fait de la variable selectionnee
        lignestemp = findall(isequal(1), A[:,bestCandidate])

        # Identifie toutes les colonnes qui doivent etre supprimées suite au candidat selectionne
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
            println("jselec : ", bestCandidate)
            println("ivar   : ", index)
            #println("C      : ", C)
            #println("A      : ", A)
            println("U      : ", candidates)
            println("-------- ")
        end
     

    end
    return sol, z

end