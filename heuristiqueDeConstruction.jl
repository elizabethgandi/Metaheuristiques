include("getfname.jl")
using (SparseArrays)

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n             = size(A,2)                   # n nombre de variables
    m             = size(A,1)                   # m nombre de contraintes
    sol           = Vector{Int64}(undef, n)     # Vecteur de base de la solution
    index         = Vector{Int64}(undef, n)     # Index d'origines des variables
    z             = 0                           # z la valeur de la fonction objective
    lignes        = Vector{Int64}(undef, n)
    colonnes      = Vector{Int64}(undef, n)
    constraints   = Vector{Vector{Int64}}(undef, m)
    constraintsLi = Vector{Vector{Int64}}(undef, n)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    for i in eachindex(constraintsLi)
        constraintsLi[i] = A[:,i]
    end

    sol        = zeros(Int, size(C, 1))
    sommeA     = zeros(Int, size(C, 1))
    index      = collect(1:size(C, 1))
    candidates = utility(C, constraints)
    lignes     = zeros(Int, size(C, 1))
    colonnes   = zeros(Int, size(C, 1))
    lignestemp = zeros(Int, size(C, 1))

    while (size(A,1) != 0) && (size(C,1) != 0) 

        # Calcul la valeur de chancune des lignes de la matrice A
        sommeA = vec(sum(A, dims=2))

        # Si il n'y a qu'une seule variable dans une des lignes elle est choisie puis supprimée
        isAlone = findfirst(x->x==minimum(sommeA), sommeA)

        # Selection de l'indice du meilleur candidat
        bestCandidate = findfirst(x->x==maximum(candidates), candidates)

        # Mise à jour de la base de la solution
        sol[index[bestCandidate]] = 1

        # Mise à jour de la valeur de la fonction objective à chaque itérations
        z = z + C[bestCandidate]

        # Suppression des lignes et colonnes dans le modele
        lignestemp = findall(isequal(1), A[:,bestCandidate])

        # Si une ligne complete est à 0 on la supprime
        isEqualZero = findfirst(isequal(0), sommeA)

        if isEqualZero !== nothing && sommeA[isEqualZero] == 0
            A = A[1:size(A,1) .!= isEqualZero,: ]
        end

        # Creation du tableau colonnetemp pour supprimer les colonnes qui sont à 1 lorsque seul le meilleur candidat doit etre à 1
        colonnetemp=(Int64)[]

        # Permet de concatener en 1 vecteur toutes les colones qui doivent étre supprimées
        for i in lignestemp
            colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
        end

        # On supprime dans les tableaux index et C les indices des colonnes supprimées
        index = index[setdiff(1:end, colonnetemp)]
        C     = C[setdiff(1:end, colonnetemp)]

        # On supprime dans la matrice A les colonnes qui doivent étre supprimées
        A = A[:, setdiff(1:end, colonnetemp)]

        # On supprime dans la matrice A la colonne du meilleur candidat
        A = A[setdiff(1:end, bestCandidate), :]
        
        # On enleve les contraintes déja satisfaite
        constraints = constraints[setdiff(1:end, lignestemp)]

        # On supprime des tableaux des index et de C celui du meilleur candidat
        index = index[setdiff(1:end, bestCandidate)]
        C     = C[setdiff(1:end, bestCandidate)]

        # On supprime les lignes de la matrice dont les contraintes sont déja satisfaites
        A = A[setdiff(1:end, lignestemp), :]
        
        # On supprime dans le tableau des colonnes à supprimer celle du meilleur candidat
        colonnetemp = setdiff(colonnetemp, bestCandidate)
       
        # On supprime dans le tableau des candidats potentiels l'indice de celui que l'on vient de prendre
        candidates = candidates[setdiff(1:end, bestCandidate)]

        # Calcul de la fonction d'utilite pour la prochaine iteration
        candidates = utility(C, constraints)
     
    end
    return sol, z
end

# Retourne la fonction utilité
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
end