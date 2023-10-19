
# Résolution SPP par algorithme glouton

function glouton(C::Vector{Int64}, A::Matrix{Int64})

    verbose::Bool = false

    # Initialisation
    n::Int64             = size(A,2)       # n nombre de variables
    m::Int64             = size(A,1)       # m nombre de contraintes

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    bestCandidate::Int64 = 0


    # 0) CALCUL DES UTILITES --------------------------------------------------

    candidates::Vector{Float64} = C ./ vec(sum(A, dims = 1))

    if verbose
        println("ivar   : ", index)
        #println("C      : ", C)
        #println("A      : ", A)
        println("U      : ", candidates)
        println(" ")
    end


    while (size(A,1) != 0) && (size(C,1) != 0) 


        # 1) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection de l'indice du meilleur candidat au regard de son utilite
        bestCandidate = argmax(candidates)
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]


        # 2) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

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

function GRASP(C::Vector{Int64}, A::Matrix{Int64})

    verbose::Bool = false

    # Initialisation
    n::Int64               = size(A,2)       # n nombre de variables
    m::Int64               = size(A,1)       # m nombre de contraintes
    alpha::Vector{Float64} = zeros(Int64,1)
    nbIterations::Int64    = 5               #stoppingRule
    limite::Float64        = 0

    index::Vector{Int64}   = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}     = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64               = 0               # z la valeur de la fonction objective

    bestCandidate::Int64   = 0
    alpha                  = rand(1)         # retourne un vecteur


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
    

    while (size(A,1) != 0) && (size(C,1) != 0) 

        # 0) MISE EN PLACE DE LA SOLUTION ALÉATOIRE ---------------------------

        limite = candidates[argmin(candidates)] + reduce(vcat, alpha)*(candidates[argmax(candidates)]-candidates[argmin(candidates)])
        RCL = findall(x-> (x >= limite), candidates)

        # 1) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection au hazard de l'indice d'un des candidats 
        bestCandidate = rand(1:size(RCL, 1))
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]


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
            println("RCL    : ", RCL)
            println("-------- ")
        end
     

    end
    return sol, z

end