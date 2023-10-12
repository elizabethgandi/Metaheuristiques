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
    #
    lignes   = Vector{Int64}(undef, n)
    colonnes = Vector{Int64}(undef, n)
    Acopy = sparse(A)
    @show A

    # Définition d'un vecteur de contraintes
    constraints = Vector{Vector{Int64}}(undef, m)
    constraintsLi = Vector{Vector{Int64}}(undef, n)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    for i in eachindex(constraintsLi)
        constraintsLi[i] = A[:,i]
    end

    sol        = zeros(Int, size(C, 1))
    #solB        = zeros(Int, size(C, 1))
    sommeA     = zeros(Int, size(C, 1))
    index      = collect(1:size(C, 1))
    #indexB      = collect(1:size(C, 1))
    candidates = utility(C, constraints)
    #
    lignes     = zeros(Int, size(C, 1))
    colonnes   = zeros(Int, size(C, 1))

    while (size(A,1) != 0) && (size(C,1) != 0) #!(isempty(candidates))

        # Calcul la valeur de chancune des lignes de la matrice A
        sommeA = vec(sum(A, dims=2))
        #sommeB = vec(sum(Acopy, dims=2))

        # Si il n'y a qu'une seule variable dans une des lignes elle est choisie puis supprimée
        isAlone = findfirst(x->x==minimum(sommeA), sommeA)
        #isAloneB = findfirst(x->x==minimum(sommeB), sommeB)

        if sommeA[isAlone] == 1 #|| sommeB[isAloneB] == 1
            A = A[1:size(A,1) .!= isAlone,: ]
            bestCandidate = isAlone
            #Acopy = Acopy[1:size(Acopy,1) .!= isAloneB,: ]
            #bestCandidateB = isAloneB
        else
            # Selection de l'indice du meilleur candidat
            bestCandidate = findfirst(x->x==maximum(candidates), candidates)
            #bestCandidateB = bestCandidate
        end
        #@show bestCandidateB
        #@show bestCandidate

        # La colonne du poids max est mise à 0
        #Acopy[:,bestCandidate] .= 0; #!pose un probleme!
        #@show Acopy
        
        # Mise à jour de la base de la solution
        sol[index[bestCandidate]] = 1
        #solB[indexB[bestCandidateB]] = 1
        #@show indexB[bestCandidateB]
        #@show index[bestCandidate]
        #@show index
        #@show indexB
        @show sol
        #@show solB

        # Mise à jour de la valeur de la fonction objective à chaque itérations
        z = z + C[bestCandidate]
        @show z

        # Suppression des lignes et colonnes dans le modele

        #lignestemp = findall(x->x==1, constraints[bestCandidate])
        #@show lignestemp

        lignestemp = findall(isequal(1), A[:,bestCandidate])
        @show lignestemp
        coltemp = findall(x->x==1, constraintsLi[lignestemp])
        @show coltemp
        

        # Si une ligne complete est à 0
        isEqualZero = findfirst(x->x==0, sommeA)

        if isEqualZero !== nothing && sommeA[isEqualZero] == 0
            A = A[1:size(A,1) .!= isEqualZero,: ]
        end

        constraints = constraints[setdiff(1:end, lignestemp)]

        #Acopy[:,lignestemp[1]] .= 0;

        #dropzeros!(Acopy)
        #@show Acopy

        coltemp = coltemp[setdiff(1:end, C)]
        @show coltemp

        index = index[setdiff(1:end, lignestemp)]
        @show index

        C = C[setdiff(1:end, lignestemp)]
        A = A[:, setdiff(1:end, lignestemp)]

        # Calcul de la fonction d'utilite pour la prochaine iteration
        candidates = utility(C, constraints)
        @show candidates

        candidates = candidates[setdiff(1:end, constraints)]
    
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