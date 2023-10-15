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
    allSup =  zeros(Int, size(C, 1))
    allSup2 =  zeros(Int, size(C, 1))

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

   # @show A
    #@show constraints
   # @show constraintsLi
    #@show candidates

    while (size(A,1) != 0) && (size(C,1) != 0) 

        # Calcul la valeur de chancune des lignes de la matrice A
        sommeA = vec(sum(A, dims=2))
       # @show sommeA

        # Si il n'y a qu'une seule variable dans une des lignes elle est choisie puis supprimée
        isAlone = findfirst(x->x==minimum(sommeA), sommeA)
        #@show isAlone


        # Selection de l'indice du meilleur candidat
        #if sommeA[isAlone] == 1 
            #A = A[1:size(A,1) .!= isAlone,: ]
         #   bestCandidate = isAlone
        #else
            bestCandidate = findfirst(x->x==maximum(candidates), candidates)
        #end
        #@show bestCandidate
        

        # Mise à jour de la base de la solution
        sol[index[bestCandidate]] = 1
        #@show sol

        #@show A
        # Mise à jour de la valeur de la fonction objective à chaque itérations
        z = z + C[bestCandidate]
        #@show z

        #@show A
        # Suppression des lignes et colonnes dans le modele
        lignestemp = findall(isequal(1), A[:,bestCandidate])
        #lignestemp = findfirst(isequal(1), A[:,bestCandidate])
        #@show lignestemp

        # Si une ligne complete est à 0
        isEqualZero = findfirst(isequal(0), sommeA)

        if isEqualZero !== nothing && sommeA[isEqualZero] == 0
            A = A[1:size(A,1) .!= isEqualZero,: ]
        end
        #@show isEqualZero

        colonnetemp=(Int64)[]
        for i in lignestemp
            colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
        end

        index = index[setdiff(1:end, colonnetemp)]

        C = C[setdiff(1:end, colonnetemp)]

        A = A[:, setdiff(1:end, colonnetemp)]
        A = A[setdiff(1:end, bestCandidate), :]
        
        constraints = constraints[setdiff(1:end, lignestemp)]
        index       = index[setdiff(1:end, bestCandidate)]

        C = C[setdiff(1:end, bestCandidate)]
        A = A[setdiff(1:end, lignestemp), :]
        
        colonnetemp = setdiff(colonnetemp,bestCandidate)
       
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