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
                        # Tester si la solution zCourant est meilleure que zBest
                        if (zCourant-C[i]-C[j]+C[k] > zMeilleur)
                            #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                            if ((findfirst(x->x>1,contraintesSaturees - A[:,i] - A[:,j] + A[:,k])) == nothing)
                                print("x")
                                #v10A
                                v1 = i
                                #v10b
                                v2 = j
                                #v01
                                v3 = k
                                #nouvelle valeur de z
                                zMeilleur = zCourant-C[i]-C[j]+C[k] 
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

            #Mise à jour des ensembles tabDesIndicesVariablesAZero et tabDesIndicesVariablesAUn

            #SUPPRIMER
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v1)]
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v2)]
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]

            #AJOUTER
            push!(tabDesIndicesVariablesAZero, v1)
            push!(tabDesIndicesVariablesAZero, v2)
            push!(tabDesIndicesVariablesAUn, v3)

            #Mise a jour du membre de droit contraintesSaturees
            contraintesSaturees -= A[:, v1]
            contraintesSaturees -= A[:, v2]
            contraintesSaturees += A[:, v3]

        end
        ameliorer = false
    end

    if (verbose)
        print("\n> 1-1 : ")   
    end
    
    ameliorer = true 
    while (ameliorer)
        ameliorer = false
        for i in tabDesIndicesVariablesAUn
            for k in tabDesIndicesVariablesAZero
                # Tester si la solution zCourant est meilleure que zBest
                if (zCourant-C[i]+C[k] > zMeilleur)
                    #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                    if ((findfirst(x->x>1,contraintesSaturees - A[:,i] + A[:,k])) == nothing) 
                        print("x")
                        #v10
                        v1 = i
                        #v01
                        v3 = k
                        #nouvelle valeur de z
                        zMeilleur = zCourant-C[i]+C[k] 
                        ameliorer = true
                    end
                end
            end 
        end
    
        if ameliorer
    
            #Mise a jour de xCourant et zCourant
            xCourant[v1] = 0
            xCourant[v3] = 1
            zCourant = zMeilleur
    
            #Mise à jour des ensembles tabDesIndicesVariablesAZero et tabDesIndicesVariablesAUn
    
            #SUPPRIMER
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v1)]
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]
    
            #AJOUTER
            push!(tabDesIndicesVariablesAZero, v1)
            push!(tabDesIndicesVariablesAUn, v3)
    
            #Mise a jour du membre de droit contraintesSaturees
            contraintesSaturees -= A[:, v1]
            contraintesSaturees += A[:, v3]
        end
        ameliorer = false
    end

    if (verbose)
        print("\n> 0-1 : ")
    end

    ameliorer = true 
    while (ameliorer)
        ameliorer = false
        for k in tabDesIndicesVariablesAZero
            # Tester si la solution zCourant est meilleure que zBest
            if (zCourant+C[k] > zMeilleur)
                #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                if ((findfirst(x->x>1,contraintesSaturees  + A[:,k])) == nothing) 
                    print("x")
                    #v01
                    v3 = k
                    #nouvelle valeur de z
                    zMeilleur = zCourant+C[k] 
                    ameliorer = true
                end
            end
        end         
    
        if ameliorer
    
            #Mise a jour de xCourant et zCourant
            xCourant[v3] = 1
            zCourant = zMeilleur
    
            #Mise à jour des ensembles tabDesIndicesVariablesAZero et tabDesIndicesVariablesAUn
    
            #SUPPRIMER
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]
    
            #AJOUTER
            push!(tabDesIndicesVariablesAUn, v3)
    
            #Mise a jour du membre de droit contraintesSaturees
            contraintesSaturees += A[:, v3]
        end
        ameliorer = false
    end

    println(" ")
    return xCourant, zMeilleur
end
