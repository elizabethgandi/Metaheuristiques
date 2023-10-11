include("getfname.jl")
using (SparseArrays)
using OrderedCollections

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n     = size(A,2)                   # n nombre de variables
    m     = size(A,1)                   # m nombre de contraintes
    sol   = Vector{Int64}(undef, n)     # Vecteur de base de la solution
    index = Vector{Int64}(undef, n)     # Index d'origines des variables
    z     = 0                           # z la valeur de la fonction objective

    lines    = OrderedSet{Int64}()
    column   = OrderedSet{Int64}()

    # Définition d'un vecteur de contraintes
    constraints = Vector{Vector{Int64}}(undef, m)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    sol        = zeros(Int, size(C, 1))
    sommeA     = zeros(Int, size(C, 1))
    index      = collect(1:size(C, 1))
    candidates = utility(C, constraints)

    while !(isempty(candidates)) 

        # Calcul la valeur de chancune des lignes de la matrice A
        sommeA = vec(sum(A, dims=2))

        # Si il n'y a qu'une seule variable dans une des lignes elle est choisie puis supprimée
        isAlone = findfirst(x->x==minimum(sommeA), sommeA)

        if sommeA[isAlone] == 1 
            A = A[1:size(A,1) .!= isAlone,: ]
            bestCandidate = isAlone
        else
            # Selection de l'indice du meilleur candidat
            bestCandidate = findfirst(x->x==maximum(candidates), candidates)
        end
        
        # Mise à jour de la base de la solution
        sol[index[bestCandidate]] = 1

        # Mise à jour de la valeur de la fonction objective à chaque itérations
        z = z + C[bestCandidate]

        # Suppression des lignes et colonnes dans le modele
        lines  = getline(constraints, bestCandidate) # identifie les éléments à 1 dans la contraintes
        column = getColumn(constraints, lines)

        # Si une ligne complete est à 0
        isEqualZero = findfirst(x->x==0, sommeA)

        if isEqualZero !== nothing && sommeA[isEqualZero] == 0
            A = A[1:size(A,1) .!= isEqualZero,: ]
        end

        sort!(column)
        sort!(lines)
        deleteat!(constraints, lines)
        
        for i in eachindex(constraints)
            deleteat!(constraints[i], column)
        end

        deleteat!(C, column)
        deleteat!(index, column)

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

# Retourne les lignes à supr en fonction de la matrice et de l'index du meilleur candidat

function getline(A, index)
    lines = OrderedSet{Int64}()

    for i in eachindex(A)
        if (A[i][index] != 0)
            push!(lines, i)
        end
    end
    return lines
end

# Retourne les colonnes à supr en fonction de la matrice et de l'ensemble des lignes à supr

function getColumn(A, lines)
    column = OrderedSet{Int64}()

    for i in lines
        for j in 1:length(A[1])
            if (A[i][j] != 0)
                push!(column, j)
            end
        end
    end
    return column
end
