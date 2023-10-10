include("getfname.jl")
using OrderedCollections

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n     = size(A,2)                   # n nombre de variables
    m     = size(A,1)                   # m nombre de contraintes
    sol   = Vector{Int64}(undef, n)     # Vecteur de base de la solution
    index = Vector{Int64}(undef, n)     # Index d'origines des variables
    z     = 0                           # z la valeur de la fonction objective

    lines  = OrderedSet{Int64}()
    column = OrderedSet{Int64}()

    # Définition d'un vecteur de contraintes
    constraints = Vector{Vector{Int64}}(undef, m)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end
    @show constraints

    sol   = zeros(Int, size(C, 1))
    sommeA   = zeros(Int, size(C, 1))
    index = collect(1:size(C, 1))

    candidates = utility(C, constraints)

    while !(isempty(candidates))

        sommeA = vec(sum(A, dims=2))
        @show sommeA

        #pré-traitrement: si il n'y a qu'une seule variable dans la ligne elle est choisie puis supprimée
        test = findfirst(x->x==minimum(sommeA), sommeA)

        if sommeA[test] == 1
            A = A[1:size(A,1) .!= test,: ]
            bestCandidate = test
        else
            # Selection de l'indice du meilleur candidat
            bestCandidate = findfirst(x->x==maximum(candidates), candidates)
        end
        
        # La colonne du poids max est mise à 0
        A[:,bestCandidate] .= 0;
        
        # Mise à jour de la base de la solution
        sol[index[bestCandidate]] = 1

        # Mise à jour de la valeur de la fonction objective à chaque itérations
        z = z + C[bestCandidate]

        # Suppression des lignes et colonnes dans le modele
        lines  = getline(constraints, bestCandidate)
        column = getColumn(constraints, lines)

        #si une ligne complete est à 0
        test2 = findfirst(x->x==0, sommeA)

        if test2 !== nothing && sommeA[test2] == 0
            A = A[1:size(A,1) .!= test2,: ]
        end

        sort!(column)
        sort!(lines)
        #@show lines
        deleteat!(constraints, lines)

        for i in eachindex(constraints)
            deleteat!(constraints[i], column)
        end

        #@show column
        deleteat!(C, column)
        #A = A[1:size(A,1) .!= :,column ]
        deleteat!(index, column)

        # Calcul de la fonction d'utilite pour la prochaine iteration
        candidates = utility(C, constraints) #PEPE faire un argmax

    end
    return sol, z
end

function utility(C, A)
    elt = Vector{Float64}(undef, length(C))

    for i in eachindex(C)
        occ = cptOcc(A, i)
        elt[i] = C[i]/occ
    end

    return elt
end

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
