# Amelioration d'une solution de SSP par algo de descente profonde
using(SparseArrays)

# Descente simple récursive
# pas de gestion des sommets déjà visiter !
function simpleDescent(C, A, x) 
    # Initialisation
    x = sparsevec(x)
    A = sparse(A)
    neighbors = Vector{SparseVector{Int64,Int64}}()

    #visited = Vector{SparseVector{Int64,Int64}}()
    #push!(visited, x)

    # Calcul voisinage N(x)
    neighbors = neighborhood_1(x)

    val = z(C, x)
    # Tant que neighbors non vide 
    while !isempty(neighbors)
        newx = pop!(neighbors)
        println("try : ", Array(newx))
        # si x' € neighbors tq z(x) < z(x')
        if val < z(C, newx)
            println(Array(newx), " may be better than x... ")
            if isAllowed(A, newx)
                println("Going deeper with ", Array(newx))
                return simpleDescent(C, A, newx)
            end
            println(Array(newx), " not allowed")
        end
        println("fail")
    end
    return Array(x)
end

function simpleDescent2(C, A, x) 
    # Initialisation
    x = sparsevec(x)
    A = sparse(A)
    neighbors = Vector{SparseVector{Int64,Int64}}()

    #visited = Vector{SparseVector{Int64,Int64}}()
    #push!(visited, x)

    # Calcul voisinage N(x)
    neighbors = neighborhood_2(A, x)

    val = z(C, x)
    # Tant que neighbors non vide 
    while !isempty(neighbors)
        newx = pop!(neighbors)
        println("try : ", Array(newx))
        # si x' € neighbors tq z(x) < z(x')
        if val < z(C, newx)
            println(Array(newx), " may be better than x... ")
            if isAllowed(A, newx)
                println("Going deeper with ", Array(newx))
                return simpleDescent(C, A, newx)
            end
            println(Array(newx), " not allowed")
        end
        println("fail")
    end
    return Array(x)
end

# Retourne la valeur de la fonction objectif
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

# i à 0, 2 elts à 1, retourne un vector de (au plus) 10 solutions possibles aléatoire
function kp(x, i)
    neighbors = Vector{SparseVector{Int64}}(undef, 0)
    index = findall(iszero, x)
    @show index
    x[i] = 0
    for j in 1:10
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


#=
# Vérifie que base est bien admissible
function isAllowed(A, x)
    index = findall(t->t==1, x)
    constraints = Vector{Int64}(undef, 0)
    # Pour chaque variable, regarder dans quelle contrainte elle apparait
    for i in index
        append!(constraints, findall(t->t==1, A[:,i]))
    end
    unique!(constraints)  # vient de decouvrir ! Super utile, a reutiliser pour le glouton !!!!!!

    # Verifie chacune des contraintes
    for i in constraints
        if sum(A[i,:]) > 1
            return false
        end
    end
    return true
end

=#


