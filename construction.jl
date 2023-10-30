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

# --------------------------------------------------------------------------- #
# Constructive greedy algorithm building a feasible solution

function GreedyConstruction(C_in, A_in)

    # Implementation style to highlight the different steps of the algorithm (could be optimized)

    C = deepcopy(C_in)
    A = deepcopy(A_in)
    verbose = false

    z::Int64 = 0  # valeur de la solution
    x = zeros(Int64, length(C)) # variables de la solution a priori toutes a zero
    posVar = collect(1:1:length(C)) # vecteur d'indices des variables
    ligneSaturee = Int64[]

    #println("\nConstruction gloutonne evaluee d'une solution admissible")

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
        j_selec = argmax(utilite)
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
# Deamon checking the feasibility of a solution given its vector x

function isAdmissible(C, A, x)

    vecSat = zeros(Int, size(A,1))
    vecUnit = ones(Int,size(A,1))
    z::Int64 = 0
    verbose = true
    var1 = findall(isequal(1), x[:])

    for j in var1
        vecSat = vecSat .+ A[:,j]
        z = z + C[j]
    end

    if findfirst(isequal(false), (vecSat .<= vecUnit)) != nothing
        println( "admissible : non")
        @assert false "detection solution non-admissible"
    end
    #println( "admissible : oui | som(x_i) = ", length(var1), " ; z = ", z)
    println( "admissible")
    return true
end
