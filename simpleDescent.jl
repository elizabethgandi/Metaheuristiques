using(SparseArrays)

#=
Tentative d'amélioration d'une solution par simple Descente
2 façons différents de générer le voisinage
simpleDescent prends a partir d'une solution met chaque variable à 1 à 0 et une variable à 0 à 1 (AddDrop)
simpleDescent2 prends la variable presente dans le plus de contraintes et la met à 0 et essaie ensuite de mettre 2 variables à 1 aléatoirement kp(1,2)
=#

# Descente simple 
# A partir d'une solution met chaque variable à 1 à 0 et une variable à 0 à 1
function simpleDescent(C, A, x) 
    # Initialisation
    x = sparsevec(x)
    A = sparse(A)
    neighbors = Vector{SparseVector{Int64,Int64}}()

    visited = Vector{SparseVector{Int64,Int64}}()
    push!(visited, x)

    # Calcul voisinage N(x)
    neighbors = neighborhood_1(x)
    # Valeur de reference de la fonction objectif
    val = z(C, x)

    # Tant qu'il reste des voisins
    while !isempty(neighbors)
        newx = pop!(neighbors)
        
        # Si le voisin est meilleur et pas déjà visité, on avance
        if !(newx in visited) && (val <= z(C, newx))
            if isAllowed(A, newx)
                x = newx
                val = z(C, x)
                neighbors = neighborhood_1(x)
            end
        end
        push!(visited, newx)
    end
    @show z(C,x)
    return Array(x)
end

# Descente simple 2
# prends la variable presente dans le plus de contraintes et la met à 0 et essaie ensuite de mettre 2 variables à 1 aléatoirement kp(1,2)
function simpleDescent2(C, A, x) 
    # Initialisation
    x = sparsevec(x)
    A = sparse(A)
    neighbors = Vector{SparseVector{Int64,Int64}}()

    visited = Vector{SparseVector{Int64,Int64}}()
    push!(visited, x)

    # Calcul voisinage N(x)
    neighbors = neighborhood_2(A, x)
    # Valeur de reference de la fonction objectif
    val = z(C, x)

    # Tant qu'il reste des voisins 
    while !isempty(neighbors)
        newx = pop!(neighbors)

        # Si le voisin est meilleur et pas déjà visité, on avance
        if !(newx in visited) && (val <= z(C, newx))
            if isAllowed(A, newx)
                x = newx
                val = z(C, x)
                neighbors = neighborhood_2(C,x)
            end
        end
        push!(visited, newx)
    end
    @show z(C,x)
    return Array(x)
end

# Retourne la valeur de la fonction objectif pour une solution x
function z(C, x::SparseVector{Int64})
    value = 0.0
    for i in findall(!iszero, x)
        value = value + C[i]
    end
    return value
end

# AddDrop pour chaque variable à 1 avec toutes les variables à 0
function neighborhood_1(x::SparseVector{Int64})
    neighbors = Vector{SparseVector{Int64}}(undef, 0)
    # Y a t il équivalence entre findall(!iszero, x) et findnz(x) ???
    for i in findall(!iszero, x)
        for j in findall(iszero, x)
            push!(neighbors, addDrop(copy(x),i,j))
        end
    end
    return neighbors
end

# Mouvement simple : i à 0, j à 1
function addDrop(x, i, j)
    x[i] = 0
    x[j] = 1
    # on supprime les zeros stockés
    dropzeros!(x)
    return x
end

# Mise à 0 de la variable apparaissant dans le plus de contraintes et mise à 1 de 2 autres contraintes nulles (k-p de type 1-2)
function neighborhood_2(A::SparseMatrixCSC{Int64, Int64}, x::SparseVector{Int64})
    coef = sum(A, dims=2)
    max = findfirst(x->x==maximum(coef), coef)
    return kp(copy(x), max)
end

# i à 0, 2 elts à 1, retourne un vecteur de (au plus) le nb de zero dans x solutions possibles aléatoire
function kp(x, i)
    neighbors = Vector{SparseVector{Int64}}(undef, 0)
    index = findall(iszero, x)
    x[i] = 0
    dropzeros!(x)
    n = max(10, length(index))
    for j in 1:n
        t = copy(x)
        t[rand(index)] = 1
        t[rand(index)] = 1
        push!(neighbors, t)
    end
    return unique!(neighbors)
end


# Vérifie que base est bien admissible
function isAllowed(A::SparseMatrixCSC{Int64, Int64}, x::SparseVector{Int64})
    index = findall(!iszero, x)
    constraints = Vector{Int64}(undef, 0)
    # Pour chaque variable, regarder dans quelles contraintes elle apparait
    for i in index
        append!(constraints, findall(t->t==1, A[:,i]))
    end
    unique!(constraints)
    # Verifie chacune des contraintes
    for i in constraints
        somme::Int64 = 0
        for j in index
            somme = somme + A[i,j]
        end
        if somme > 1
            return false
        end
    end
    return true
end
