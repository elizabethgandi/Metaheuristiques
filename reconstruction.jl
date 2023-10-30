# --------------------------------------------------------------------------- #
# Reconstruction façon VND

function reconstructionDR(C, A, x_in, z_in)

    verbose = false

    xBest::Vector{Int}   = copy(x_in) # meilleure solution
    zBest::Int           = z_in
    zCourant::Int        = z_in       # solution courante

    # Split le vecteur x[i] en deux vecteurs, contenant les indices des x[i] à 0 et x[i] à 1
    var0::Vector{Int}    = findall(isequal(0), x_in[:]) # variables a zero
    var1::Vector{Int}    = findall(isequal(1), x_in[:]) # variables a un

    # Resume dans contraintesSaturees les contraintes qui sont saturees par la solution x[i]
    contraintesSaturees::Vector{Int} = zeros(Int, size(A,1)) 
    for j in var1
        contraintesSaturees = contraintesSaturees + A[:,j] 
    end

    nL = size(A,1)

    # 0-1 exchange ------------------------------------------------------------

    print("      > 0-1 : ")

    ameliorer = true
    while ameliorer == true
        
        # recherche plus profonde descente dans le voisinage de la solution courante avec un mouvement 0-1 exchange
        j01Best = 0
        ameliorer = false
        for j01 in var0
            # test si ameliore zbest actuel
            # test l'admissibilite de l'ajout d'une variable a la solution
            
            if  (zCourant + C[j01] > zBest) 
                valide = true
                i = 1
                while valide && i≤nL
                    if contraintesSaturees[i] + A[i,j01] > 1
                        valide = false
                    else
                        i+=1
                    end
                end
                if valide    
                    # meilleur et aucune contrainte violee => meilleur voisin admissible
                    j01Best = j01
                    zBest = zCourant + C[j01Best]
                    ameliorer = true
                end
            end          
        end
        if ameliorer
            # maj meilleure solution courante
            var0 = setdiff(var0,j01Best)     # maj var0
            push!(var1,j01Best)              # maj var1
            xBest[j01Best] = 1               # maj xBest
            zCourant = zBest                 # maj zCourant
            contraintesSaturees = contraintesSaturees + A[:,j01Best] # maj membre de droite

            print("x")
        end
    end
    println(" ")

    # 1-1 exchange ------------------------------------------------------------

    print("      > 1-1 : ")

    ameliorer = true
    while ameliorer == true

        # recherche plus profonde descente dans le voisinage de la solution courante avec un mouvement 1-1 exchange

        j01Best = 0
        j10Best = 0
        ameliorer = false
    
        for j10 in var1
            for j01 in var0
                # test si ameliore zbest actuel
                # test l'admissibilite de l'ajout d'une variable a la solution
    
                if  (zCourant - C[j10] + C[j01] > zBest) 
                    valide = true
                    i = 1
                    while valide && i≤nL
                        if contraintesSaturees[i] - A[i,j10] + A[i,j01] > 1
                            valide = false
                        else
                            i+=1
                        end
                    end
                    if valide
                        # meilleur et aucune contrainte violee => meilleur voisin admissible
                        j01Best = j01
                        j10Best = j10
                        zBest = zCourant - C[j10Best] + C[j01Best]
                        ameliorer = true
                    end                    
                end
            end
        end
        if ameliorer
            # maj meilleure solution courante
            var0 = setdiff(var0,j01Best)                    # maj var0
            push!(var0,j10Best)           
            var1 = setdiff(var1,j10Best)                    # maj var1
            push!(var1,j01Best)           
            xBest[j01Best] = 1                              # maj xBest
            xBest[j10Best] = 0
            zCourant = zBest                                # maj zCourant
            contraintesSaturees = contraintesSaturees - A[:,j10Best] + A[:,j01Best] # maj membre de droite

            print("x")
        end
    end
    println("")

    # 2-1 exchange ------------------------------------------------------------

    print("      > 2-1 : ")

    ameliorer = true
    while ameliorer == true
        # recherche plus profonde descente dans le voisinage de la solution courante avec un mouvement 2-1 exchange

        j01Best  = 0
        j10aBest = 0
        j10bBest = 0
        var1a    = copy(var1)
        var1b    = copy(var1)
        ameliorer   = false

        for j10a in var1a
            for j10b in setdiff(var1b,j10a)
                for j01 in var0
                    # test si ameliore zbest actuel
                    # test l'admissibilite de l'ajout d'une variable a la solution

                    if  (zCourant - C[j10a] - C[j10b] + C[j01] > zBest) 
                        valide = true
                        i = 1
                        while valide && i≤nL
                            if contraintesSaturees[i] - A[i,j10a] - A[i,j10b] + A[i,j01] > 1
                                valide = false
                            else
                                i+=1
                            end
                        end
                        if valide
                            # meilleur et aucune contrainte violee => meilleur voisin admissible
                            j01Best = j01
                            j10aBest = j10a
                            j10bBest = j10b
                            zBest = zCourant - C[j10aBest] - C[j10bBest] + C[j01Best]
                            ameliorer = true
                        end
                    end
                end
            end
        end
        if ameliorer
            # maj meilleure solution courante
            var0 = setdiff(var0,j01Best)                                     # maj var0
            push!(var0,j10aBest)
            push!(var0,j10bBest)
            var1 = setdiff(var1,j10aBest)                                    # maj var1
            var1 = setdiff(var1,j10bBest)
            push!(var1,j01Best)
            xBest[j01Best] = 1                                               # maj xBest
            xBest[j10aBest] = 0
            xBest[j10bBest] = 0
            zCourant = zBest                                                 # maj zCourant
            contraintesSaturees = contraintesSaturees - A[:,j10aBest] - A[:,j10bBest] + A[:,j01Best] # maj membre de droite

            print("x")
        end
    end

    println("")

    return xBest, zBest
end
