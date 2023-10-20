
include("getfname.jl")
# Ne plus utiliser les sparses Array
# Evaluer le cout de la fonction avant même de créer le voisin
# Evaluer l'admissibilité de la fonction avant même de créer le voisin sans tout regarder
# Créer le voisin meilleur et recommencer la recherche sur lui avec le même Mouvement
# Passer au mouvement suivant

# Créer un vecteur contenant les indices de 0 pour pouvoir itérer dessus
# Même chose pour les 1

# 3 voisinages
# 1:2
# 1:1
# 0:1


# Descente simple
# A partir d'une solution met chaque variable à 1 à 0 et une variable à 0 à 1
function simpleDescent(C::Vector{Int64}, A::Matrix{Int64}, x::Vector{Int64}, z::Int64) 

    # Affichage 
    verbose::Bool = true

    # Init
    bestX               = zeros(Int64, length(C))
    bestZ               = z
    neighbor            = x
    admissible::Bool    = true
    l::Int64            = 1

    var0::Vector{Int64} = findall(iszero, x)                # Indice des variables à 0
    var1::Vector{Int64} = findall(!iszero, x)               # Indice des variables à 1
    ctr::Vector{Int64}  = zeros(Int64, size(A, 1))          # Contraintes utilisées

    for i in var1
        for j in findall(!iszero, A[:,i])
            ctr[j] = 1
        end
    end

    #println("Contraintes actives :", ctr)

    # Boucle mouvement 2:1
    for i in var0
        # mettre i à 1
        for j in var0
            # mettre j à 1
            for z in var1
                #mettre z à 0
                # verif meilleur
                if bestZ < (bestZ - C[z] + C[j] + C[i])
                    #println("Vérification d'admissibilité !")

                    l = 1
                    while ((l < size(A,1)) && (admissible))
                    #for l in 1:size(A,1)
                        # + valeur de contrainte active 
                        # + somme des valeurs des variables mise à 1 dans cette contrainte
                        # - somme des valeurs des variables mise à 0 dans cette contrainte
                        # si > 1 alors + de deux variables sont actives sur la contrainte, solution non admissible
                        if (ctr[l] + A[l,i] + A[l,j] - A[l,z]) > 1
                            admissible = false
                        end
                        l = l+1
                    end

                    if admissible
                    # vérif d'admissibilité
                    #if isChangeAdmissible2!(A, [i,j], [z], ctr)
                        println("Meilleur trouvé !")
                        # Mise à jour de bestZ, var0, var1
                        bestZ = bestZ + C[i] + C[j] - C[z]
                        append!(var0, z)
                        deleteat!(var0, findall(m->(m==i||m==j), var0)) #Remplace par un setdiff ?
                        append!(var1, [i,j])
                        deleteat!(var1, findall(m->m==z, var1))

                        #update ctr

                        # Conraintes à mettre à 0
                        update::Vector{Int64} = findall(!iszero, A[:,z])
                        for m in update
                            ctr[m] = 0
                        end

                        update = findall(!iszero, A[:,z])
                        for m in update
                            ctr[m] = 1
                        end

                        update = findall(!iszero, A[:,z])
                        for m in update
                            ctr[m] = 1
                        end

                        # joli affichage
                    end
                    admissible = true
                end 
            end    
        end
    end

    # Construction de la solution
    for i in var1
        bestX[i] = 1
    end

    return bestX, bestZ
end


# Fonction "générale" pour un k-p exchange
# Suppose que la solution était déjà admissible avant modification
# A matrice des contraintes
# to1 vector de variable mise à 1
# to0 vector de variable mise à 0
# var1 vector de variable déjà à 1

function isChangeAdmissible!(A::Matrix{Int64}, to1::Vector{Int64}, to0::Vector{Int64}, ctr::Vector{Int64})
    # Init
    # Vérifie que les contraintes sont respectées
    for i in 1:size(A,1)
        # + valeur de contrainte active (1 si contrainte déjà activé par une variable 0 sinon)
        # + somme des valeurs des variables mise à 1 dans cette contrainte (1 si active 0 sinon pour chaque variable)
        # - somme des valeurs des variables mise à 0 dans cette contrainte (1 si active 0 sinon pour chaque variable)
        # si > 1 alors + de deux variables sont actives sur la contrainte, solution non admissible
        if (ctr[i] + sum(A[i,to1[k]] for k in 1:length(to1)) - sum(A[i,to0[k]] for k in 1:length(to0))) > 1
            return false
        end
    end
    return true
end