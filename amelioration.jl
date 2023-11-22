function gloutonAmelioration(C, A, xConstruction, zConstruction)

    # Intinialisation -----------------------------------------------

    xCourant                    = copy(xConstruction)
    tabDesIndicesVariablesAZero = findall(isequal(0), xConstruction) # tableau des indices des variables qui sont à 0 dans xConstruction
    tabDesIndicesVariablesAUn   = findall(isequal(1), xConstruction) # tableau des indices des variables qui sont à 1 dans xConstruction

    zCourant                    = zConstruction
    zMeilleur                   = zConstruction
    v1                          = 0                                  # indice i la première variable mise à 0, si la solution z est meilleure
    v2                          = 0                                  # indice j la seconde variable mise à 0, si la solution z est meilleure
    v3                          = 0                                  # indice k la variable mise à 1, si la solution z est meilleure

    ameliorer                   = true                               # permet de savoir si notre solution trouvée est améliorante ou non, auquel cas on continue dans la boucle tant que
    verbose                     = false                              # utilisé pour les affichages, faux -> aucun affichages, vrai -> tous les affichages

    # 1) CALCUL DU MEMBRE DE DROITE DES CONTRAINTES -----------------
    contraintesSaturees         = zeros(Int64, size(A,1)) 
    for i in tabDesIndicesVariablesAUn
        contraintesSaturees += A[:, i]
    end
    
    # I/2) 2-1 ECHANGE ------------------------------------------------

    if (verbose)
        print("\n> 2-1 : ")
    end

    while (ameliorer)
        ameliorer = false
        for i in tabDesIndicesVariablesAUn
            for j in tabDesIndicesVariablesAUn 
                if (i<j)
                    for k in tabDesIndicesVariablesAZero

                        # 3) TESTER SI LA SOLUTION zCourant EST MEILLEURE QUE zBest
                        if (zCourant-C[i]-C[j]+C[k] > zMeilleur)

                            # 4) REMPLISSAGE DU MEMBRE DE DROITE ET VERIFICATION D'ADMISSIBILITE
                            if ((findfirst(x->x>1,contraintesSaturees - A[:,i] - A[:,j] + A[:,k])) == nothing)
                                if (verbose) 
                                    print("x") 
                                end
                                # v10A: première variable (A) qui était à 1 et qui est mise à 0 dans notre solution 
                                v1 = i
                                # v10b: seconde variable  (B) qui était à 1 et qui est mise à 0 dans notre solution 
                                v2 = j
                                # v01: unique variable qui était à 0 et qui est mise à 1 dans notre solution 
                                v3 = k

                                # Nouvelle valeur de zMeilleur
                                zMeilleur = zCourant-C[i]-C[j]+C[k] 

                                # On a améliorer notre solution alors ameliorer = true
                                ameliorer = true
                            end
                        end
                    end 
                end
            end
        end

        if ameliorer

            # 5) MISE A JOUR DE xCourant ET zCourant
            xCourant[v1] = 0
            xCourant[v2] = 0
            xCourant[v3] = 1
            zCourant = zMeilleur

            # 6) MISE A JOUR DES ENSEMBLES tabDesIndicesVariablesAZero ET tabDesIndicesVariablesAUn

            # 6.1) SUPPRIMER
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v1)]
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v2)]
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]

            # 6.2) AJOUTER
            push!(tabDesIndicesVariablesAZero, v1)
            push!(tabDesIndicesVariablesAZero, v2)
            push!(tabDesIndicesVariablesAUn, v3)

            # 7) MISE A JOUR DU MEMBRE DE DROITE contraintesSaturees
            contraintesSaturees -= A[:, v1]
            contraintesSaturees -= A[:, v2]
            contraintesSaturees += A[:, v3]

        end
    end


    # II/2) 1-1 ECHANGE ------------------------------------------------

    if (verbose)
        print("\n> 1-1 : ")   
    end
    
    ameliorer = true 
    while (ameliorer)
        ameliorer = false
        for i in tabDesIndicesVariablesAUn
            for k in tabDesIndicesVariablesAZero

                # 3) TESTER SI LA SOLUTION zCourant EST MEILLEURE QUE zBest
                if (zCourant-C[i]+C[k] > zMeilleur)
                    
                    # 4) REMPLISSAGE DU MEMBRE DE DROITE ET VERIFICATION D'ADMISSIBILITE
                    if ((findfirst(x->x>1,contraintesSaturees - A[:,i] + A[:,k])) == nothing) 
                        if (verbose) 
                            print("x") 
                        end
                        
                        # v10A: unique variable qui était à 1 et qui est mise à 0 dans notre solution 
                        v1 = i
                        # v01: unique variable qui était à 0 et qui est mise à 1 dans notre solution 
                        v3 = k
                        
                        # Nouvelle valeur de zMeilleur
                        zMeilleur = zCourant-C[i]+C[k] 
                        
                        # On a améliorer notre solution alors ameliorer = true
                        ameliorer = true
                    end
                end
            end 
        end
    
        if ameliorer
    
            # 5) MISE A JOUR DE xCourant ET zCourant
            xCourant[v1] = 0
            xCourant[v3] = 1
            zCourant = zMeilleur
    
            # 6) MISE A JOUR DES ENSEMBLES tabDesIndicesVariablesAZero ET tabDesIndicesVariablesAUn

            # 6.1) SUPPRIMER
            tabDesIndicesVariablesAUn = tabDesIndicesVariablesAUn[setdiff(1:end, v1)]
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]
    
            # 6.2) AJOUTER
            push!(tabDesIndicesVariablesAZero, v1)
            push!(tabDesIndicesVariablesAUn, v3)
    
            # 7) MISE A JOUR DU MEMBRE DE DROITE contraintesSaturees
            contraintesSaturees -= A[:, v1]
            contraintesSaturees += A[:, v3]
        end
    end

    # III/2) 0-1 ECHANGE ------------------------------------------------

    if (verbose)
        print("\n> 0-1 : ")
    end

    ameliorer = true 
    while (ameliorer)
        ameliorer = false
        for k in tabDesIndicesVariablesAZero
           
            # 3) TESTER SI LA SOLUTION zCourant EST MEILLEURE QUE zBest
            if (zCourant+C[k] > zMeilleur)
                
                # 4) REMPLISSAGE DU MEMBRE DE DROITE ET VERIFICATION D'ADMISSIBILITE
                if ((findfirst(x->x>1,contraintesSaturees  + A[:,k])) == nothing) 
                    if (verbose) 
                        print("x") 
                    end
                    
                    # v01: unique variable qui était à 0 et qui est mise à 1 dans notre solution 
                    v3 = k
                    
                    # Nouvelle valeur de zMeilleur
                    zMeilleur = zCourant+C[k] 
                    
                    # On a améliorer notre solution alors ameliorer = true
                    ameliorer = true
                end
            end
        end         
    
        if ameliorer
    
            # 5) MISE A JOUR DE xCourant ET zCourant
            xCourant[v3] = 1
            zCourant = zMeilleur
    
            # 6) MISE A JOUR DE L' ENSEMBLE tabDesIndicesVariablesAZero

            # 6.1) SUPPRIMER
            tabDesIndicesVariablesAZero = tabDesIndicesVariablesAZero[setdiff(1:end, v3)]
    
            # 6.2) AJOUTER
            push!(tabDesIndicesVariablesAUn, v3)
    
            # 7) MISE A JOUR DU MEMBRE DE DROITE contraintesSaturees
            contraintesSaturees += A[:, v3]
        end
    end

    println(" ")

    return xCourant, zMeilleur
end
