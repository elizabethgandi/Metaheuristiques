include("getfname.jl")

# Résoudre pb de division par 0

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n = size(A,2)       # n nombre de variables
    @show(n)
    m = size(A,1)       # m nombre de contraintes
    sol = Vector{Int64}(undef, n)     # Vecteur de base de la solution
    index = Vector{Int64}(undef, n)

    # Définition d'un vecteur de contraintes
    constraints = Vector{Vector{Int64}}(undef, m)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    for i in eachindex(index)
        index[i] = i
        sol[i] = 0
    end
    @show(index)

    candidates = utility(C, constraints)

    while cond(candidates)
        # Estimation du meilleur candidat avec la fonction d'utilite
        
        @show candidates
        # Selection de l'indice du meilleur candidat

        bestCandidate = findfirst(x->x==maximum(candidates), candidates)
        #Mauvais indice !!

        @show(bestCandidate)
        @show(index)
        sol[index[bestCandidate]] = 1

        # Mise a jour des coefficents et contraintes (? mise a 0 ?)
        # Coefficients

        deleteat!(C, bestCandidate)
        deleteat!(index, bestCandidate)
        # Contraintes
        for i in 1:length(constraints)
            deleteat!(constraints[i], bestCandidate)
        end
        candidates = utility(C, constraints)
    end
    @show sol
    return sol
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

function cond(v)
    for i in eachindex(v)
        if v[i] > 0
            return true
        end
    end
    return false
end
