function gloutonAmelioration(C, A, xConstruction, zConstruction)
    xCourant = copy(xConstruction)
    tabDesIndicesVariablesAZero = findall(isequal(0), xConstruction)
    tabDesIndicesVariablesAUn   = findall(isequal(1), xConstruction)

    zCourant  = zConstruction
    zMeilleur = zConstruction
    v1 = 0
    v2 = 0
    v3                  = 0
    ameliorer           = true
    verbose             = true

    contraintesSaturees = zeros(Int64, size(A,1)) #RHS

    for i in tabDesIndicesVariablesAUn
        contraintesSaturees += A[:, i]
    end
    #@show contraintesSaturees
    
    # 2-1 echange --------------------------------------------------

    if (verbose)
        print("\n> 2-1 : ")
    end

    while (ameliorer)
        ameliorer = false
        for i in tabDesIndicesVariablesAUn
            for j in tabDesIndicesVariablesAUn 
                if (i<j)
                    for k in tabDesIndicesVariablesAZero
                        #println(" i = $i j = $j k = $k")

                        # Tester si la solution zCourant est meilleure que zBest
                        if (zCourant-C[i]-C[j]+C[k] > zMeilleur)
                            #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                            if ((findfirst(x->x>1,contraintesSaturees - A[:,i] - A[:,j] + A[:,k])) == nothing) #&& isAdmissible(xAmelioration) == true)
                                println("x")
                                println("un voisin ameliorant trouve")
                                #v10A
                                v1 = i
                                #v10b
                                v2 = j
                                #v01
                                v3 = k
                                #nouvelle valeur de z
                                zMeilleur = zCourant-C[i]-C[j]+C[k] 
                                println("v1=$v1 v2=$v2 v3=$v3 zmeil=$zMeilleur")
                                ameliorer = true
                            end
                        end
                    end 
                end
            end
        end
        if ameliorer
            #Mise a jour de xCourant et zCourant
            xCourant[v1] = 0
            xCourant[v2] = 0
            xCourant[v3] = 1
            zCourant = zMeilleur
            println("xCourant $xCourant et zCourant $zCourant")

            #Mise à jour des ensembles tabDesIndicesVariablesAZero et tabDesIndicesVariablesAUn

            #SUPPRIMER
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v1)]
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v2)]

            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]

            #println("LORS DE L'AJOUT: var1 $tabDesIndicesVariablesAUn et var0 $tabDesIndicesVariablesAZero")

            #AJOUTER
            push!(tabDesIndicesVariablesAZero, v1)
            push!(tabDesIndicesVariablesAZero, v2)

            push!(tabDesIndicesVariablesAUn, v3)

            #println("LORS DE LA SUPPRESSION: var1 $tabDesIndicesVariablesAUn et var0 $tabDesIndicesVariablesAZero")

            #Mise a jour du membre de droit contraintesSaturees
            contraintesSaturees -= A[:, v1]
            contraintesSaturees -= A[:, v2]

            contraintesSaturees += A[:, v3]

           #println("RHS $contraintesSaturees")
        end

        ameliorer = false
    end

    @show v1, v2, v3, zMeilleur
    #@show C
    #@show xConstruction
    #println(A)

end

#=
    if (verbose)
        print("\n> 1-1 : ")
    end

    while (ameliorer)
        ameliorer = false
        for i in tabDesIndicesVariablesAUn
             for k in tabDesIndicesVariablesAZero
                        # Tester si la solution zCourant est meilleure que zBest
                        if (zCourant-C[i]+C[k] > zMeilleur)
                            #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                            if ((findfirst(x->x>1,contraintesSaturees - A[:,i] + A[:,k])) == nothing) #&& isAdmissible(xAmelioration) == true)
                                println("x")
                                println("un voisin ameliorant trouve dans 1-1")
                                #v10A
                                v1 = i
                                #v01
                                v3 = k
                                #nouvelle valeur de z
                                zMeilleur = zCourant-C[i]+C[k] 
                                println("v1=$v1 v3=$v3 zmeil=$zMeilleur")
                                #ameliorer = true
                            end
                        end
                    end 
                end
            end
        end
    end 

     if (verbose)
        print("\n> 0-1 : ")
    end

    while (ameliorer)
        ameliorer = false
                for k in tabDesIndicesVariablesAZero
                    # Tester si la solution zCourant est meilleure que zBest
                    if (zCourant+C[k] > zMeilleur)
                        #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                        if ((findfirst(x->x>1,contraintesSaturees + A[:,k])) == nothing) #&& isAdmissible(xAmelioration) == true)
                            println("x")
                            println("un voisin ameliorant trouve dans 0-1")
                            #v01
                            v3 = k
                            #nouvelle valeur de z
                            zMeilleur = zCourant+C[k] 
                            println("v1=$v1 zmeil=$zMeilleur")
                            #ameliorer = true
                        end
                    end
                end 
            end
        end
    end 
=#


#function isAdmissible(xSolution)
#
#end

function isAd(C, A, x)

    vecSat = zeros(Int, size(A,1))
    vecUnit = ones(Int,size(A,1))
    z::Int64 = 0
    verbose = true
    var1 = findall(isequal(1), x[:])
    
    for j in var1
        vecSat = vecSat .+ A[:,j]
        z = z + C[j]
    end
    
    if findfirst(isequal(false), (vecSat .<= vecUnit)) != nothing
        println( "admissible : non")
        @assert false "detection solution non-admissible"
    end
    println( "admissible : oui | som(x_i) = ", length(var1), " ; z = ", z)
    return true

end

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