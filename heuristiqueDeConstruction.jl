include("getfname.jl")
using (SparseArrays)

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n             = size(A,2)                   # n nombre de variables
    m             = size(A,1)                   # m nombre de contraintes
    #sol           = Vector{Int64}(undef, n)     # Vecteur de base de la solution
    #index         = Vector{Int64}(undef, n)     # Index d'origines des variables
    z             = 0                           # z la valeur de la fonction objective
    #lignes        = Vector{Int64}(undef, n)
    #colonnes      = Vector{Int64}(undef, n)
    constraints   = Vector{Vector{Int64}}(undef, m)
    constraintsLi = Vector{Vector{Int64}}(undef, n)
    verbose = true

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    for i in eachindex(constraintsLi)
        constraintsLi[i] = A[:,i]
    end

    sol        = zeros(Int, size(C, 1))
    sommeA     = zeros(Int, size(C, 1))
    index      = collect(1:size(C, 1))

    #candidates = utility(C, constraints)
    #lignes     = zeros(Int, size(C, 1))
    #colonnes   = zeros(Int, size(C, 1))
    lignestemp = zeros(Int, size(C, 1))
    #colonnetemp = zeros(Int, size(C, 1))
    bestCandidate = 0

    # calcul des utilites
    candidates = C ./ vec(sum(A, dims = 1))
    @show candidates

    #for boucle = 1:2 #
        while (size(A,1) != 0) && (size(C,1) != 0) 

        # Calcul la somme de chacune des lignes de A
        sommeA = vec(sum(A, dims=2))

        # INUTILE Si il n'y a qu'une seule variable dans une des lignes elle est choisie puis supprimée
        #isAlone = findfirst(x->x==minimum(sommeA), sommeA)


        # 1) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection de l'indice du meilleur candidat au regard de son utilite
        bestCandidate = argmax(candidates)
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]

        if verbose
            #@show boucle
            println("ivar   : ", index)
            #println("C      : ", C)
            #println("A      : ", A)
            println("U      : ", candidates)
            println("jselec : ",bestCandidate)
            println("-------- ")
        end

        # 2) REDUCTION DU PROBLEME SUIOTE AU CANDIDAT SELECTIONNE -------------

        # Identification des contraintes a supprimer du fait de la variable selectionnee
        lignestemp = findall(isequal(1), A[:,bestCandidate])

        # Si une ligne complete est à 0 on la supprime
        #isEqualZero = findfirst(isequal(0), sommeA)

        #if isEqualZero !== nothing && sommeA[isEqualZero] == 0
        #    A = A[1:size(A,1) .!= isEqualZero,: ]
        #end


        # identifie toutes les colonnes qui doivent etre supprimées suite au candidat selectionne
        colonnetemp=(Int64)[]
        for i in lignestemp
            # scrute contrainte par contrainte les coefficients de A de valeur 1 (et elimine les eventuels doublons)
            colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
        end

        # On supprime dans les structures index, A et C les valeurs corrspondantes aux colonnes supprimées
        index      = index[setdiff(1:end, colonnetemp)]      # indice de variables
        A          = A[:, setdiff(1:end, colonnetemp)]       # matrice des contraines
        C          = C[setdiff(1:end, colonnetemp)]          # coefficients de la fonction objectif
        candidates = candidates[setdiff(1:end, colonnetemp)] #

        # On supprime dans la matrice A les lignes correspondantes aux contraintes impliquees par la variable selectionnee
        A      = A[setdiff(1:end, lignestemp), :]
        sommeA = sommeA[setdiff(1:end, lignestemp)]

        # INUTILE ??? On enleve les contraintes déja satisfaite
        #constraints = constraints[setdiff(1:end, lignestemp)]

        # INUTILE??? On supprime des tableaux des index et de C celui du meilleur candidat
        #index = index[setdiff(1:end, bestCandidate)]
        #C     = C[setdiff(1:end, bestCandidate)]

        # INUTILE ??? On supprime les lignes de la matrice dont les contraintes sont déja satisfaites
        #A = A[setdiff(1:end, lignestemp), :]
        
        # INUTILE ??? On supprime dans le tableau des colonnes à supprimer celle du meilleur candidat
        #colonnetemp = setdiff(colonnetemp, bestCandidate)
       
        # On supprime dans le tableau des candidats potentiels l'indice de celui que l'on vient de prendre
        #candidates = candidates[setdiff(1:end, bestCandidate)]


       #= if verbose
           println("ivar   : ", colonnetemp)
            #println("C      : ", C)
            #println("A      : ", A)
            println("U      : ", candidates)
            println("jselec : ",bestCandidate)
            println("-------- ")
        end=#


        # Calcul de la fonction d'utilite pour la prochaine iteration
       # candidates = utility(C, constraints)

        #@assert stop
     
    end
    return sol, z
end

#= Retourne la fonction utilité
function utility(C, A)
    elt = Vector{Float64}(undef, length(C))

    for i in eachindex(C)
        occ = cptOcc(A, i)
        elt[i] = C[i]/occ
    end

    return elt
end

# Retourne l' occurence de chacune des variables 
function cptOcc(v::Vector{Vector{Int64}}, column)
    cpt = 0

    for i in eachindex(v)
        cpt = cpt + v[i][column]
    end

    return cpt
end=#