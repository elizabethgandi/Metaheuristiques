
include("getfname.jl")


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

    # Initialiation du vecteur de contraintes actives
    for i in var1
        for j in findall(!iszero, A[:,i])
            ctr[j] = 1
        end
    end

    if verbose
        println("N with 2:1 exchange (# number of upgrade):")
    end

    # Voisinage avec mouvement de 2:1 exchange
    for i in var0
        for j in var0
            for z in var1
                if j>i  # permet d'éviter les doublons (i,j) et (j,i) et les (i,i)
                    if bestZ < (bestZ - C[z] + C[j] + C[i])

                        # Vérification d'admissibilié de la solution
                        l = 1
                        while ((l <= size(A,1)) && (admissible))
                            if (ctr[l] + A[l,i] + A[l,j] - A[l,z]) > 1
                                admissible = false
                            end
                            l = l+1
                        end

                        if admissible # si admissible

                            if verbose 
                                print("#")
                            end

                            # Mise à jour de bestZ, var0, var1
                            bestZ = bestZ + C[i] + C[j] - C[z]
                            append!(var0, z)
                            deleteat!(var0, findall(m->(m==i||m==j), var0))
                            append!(var1, [i,j])
                            deleteat!(var1, findall(m->m==z, var1))

                            #update ctr

                            # Contraintes occupée par z sont libérées
                            update = findall(!iszero, A[:,z])
                            for m in update
                                ctr[m] = 0
                            end

                            # Contraintes occupée par i sont ajoutées
                            update = findall(!iszero, A[:,i])
                            for m in update
                                ctr[m] = 1
                            end

                            # Contraintes occupée par j sont ajoutées
                            update = findall(!iszero, A[:,j])
                            for m in update
                                ctr[m] = 1
                            end
                        end

                        admissible = true
                    end
                end 
            end    
        end
    end

    # Solution améliorée
    for i in var1
        bestX[i] = 1
    end

    print("\ntime :")
    return bestX, bestZ
end