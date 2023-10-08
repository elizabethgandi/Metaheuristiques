# Amelioration d'une solution de SSP par algo de descente profonde

function deepestDescent(C, A, x) 

    # Calcul voisinage N(x)
    # Pour tout x' € N(x), noter que x' est visité

    # Tant qu'il existe un x"€ N(x) tq x" meilleur ou équivalent à x et x' € N(x) et n'est pas déjà visité
        # x = x'
        # Calcul voisinage x
        # Pour tout x' € N(x), noter que x' est visité
    # end
    # retourner x
end


# Retourne la valeur de la fonction objectif
function z(C, x)
    value::Float64 = 0
    for i in eachindex(x)
        if x[i] == 1
            value = value + C[i]
        end
    end
    return value
end

# Mouvement x -> x'
#=
Idee 
met xi avec le plus de variables dans ses contraintes à 0  (k=1)
met un maximum de xj libéré à 1 (p=|xj|)
vérifier que la solution est bien admissible
=#
function mvt(x)

end

# Vérifie que base est bien admissible
function isAllowed(x, A)
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

