include("getfname.jl")

# Résoudre pb de division par 0

# Résolution SPP par algorithme glouton
function glouton(C, A)

    # Initialisation
    n = size(A,2)       # n nombre de variables
    m = size(A,1)       # m nombre de contraintes
    sol = fill(0, length(C));     # Vecteur de base de la solution
    
    # Définition d'un vecteur de contraintes
    constraints = Vector{Vector{Int64}}(undef, m)

    for i in eachindex(constraints)
        constraints[i] = A[i,:]
    end

    candidates = utility(C, constraints)

    while cond(candidates)
        # Estimation du meilleur candidat avec la fonction d'utilite
        
        @show candidates
        # Selection de l'indice du meilleur candidat

        bestCandidate = findfirst(x->x==maximum(candidates), candidates)
        #Mauvais indice !!


        @show(bestCandidate)
        sol[bestCandidate] = 1

        # Mise a jour des coefficents et contraintes (? mise a 0 ?)
        # Coefficients
        #C[bestCandidate] = 0
        deleteat!(C, bestCandidate)
        # Contraintes
        for i in 1:length(constraints)
            deleteat!(constraints[i], bestCandidate)
            #=
            if A[i, bestCandidate] != 0
                for j in 1:size(A,2)
                    A[i,j] = 0
                end
                deleteat!(A,i)
            end
            =#
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
