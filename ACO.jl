function ACO(C, A, nbIterationsACO, nbFourmis)
    # Initialisation

    #print(C)

    vecteurSol::Vector{Int64}          = zeros(length(C))
    vecteurPropa::Vector{Int64}        = zeros(length(C))
    vecteurPheromones::Vector{Float64} = Vector(undef, length(C)) #proba entre 0 et 1

    cheminFourmis::Vector{Vector{Int64}}                = fill([], nbFourmis) #length(C))
    vecteurSolFourmis::Vector{Int64}                    = zeros(nbFourmis)
    vecteurDeToutesLesSolutionsFinales::Vector{Float64} = zeros(length(C))

    idbestSol::Int64  = 0 
    lambda::Float64   = 0.1  # coef evaporation des fourmis initié à 10%
    beta::Float64     = 1.1  # coef evaporation des fourmis
    meilleur::Float64 = 0
    k::Int64          = 0
    z::Int64          = 0
    resFinal::Int64   = 0

    for i in 1:length(C)
        vecteurPheromones[i] = 1/length(C)
    end

    # Premier appel pour obtenir une fonction d'utilité des phéromones
    #vecteurSol = gloutonConstruction(C, A)

    # Boucle POUR pour avoir n itérations
    for i in 1:nbIterationsACO # fixer nbIterationsACO à 1 pour l'instant?

       # println("LANCÉ $i --------------------------------------------------")

        idbestSol = 0
        k = k+1

        # Boucle POUR utilisée pour chaque fourmis
        for j in 1:nbFourmis
            #println("FOURMIS $j --------------------------------------------------")
            # Fonction qui construit une solution avec la roulette
            cheminFourmis[j], z = fourmisConstruction(vecteurSol, vecteurPheromones, C, A)
            #@show cheminFourmis
            
            vecteurSolFourmis[j] = sum(cheminFourmis[j][m]*C[m] for m in 1:length(C))
            #println("vect sol $vecteurSolFourmis")
            #println("-------------------------------------------------------------")
        end

        idbestSol = argmax(vecteurSolFourmis) 

        # jusque ici tout ok-------------------------------

        # Evaporation
        for j in 1:length(C)
            #if (j == idbestSol)
                #vecteurPheromones[j] = (vecteurPheromones[j] * lambda)*5
            #else 
                vecteurPheromones[j] = vecteurPheromones[j] * lambda
           # end     
        end

        #@show vecteurPheromones

        #vecteurPheromones = evaporation(vecteurPheromones, lambda)??? Pck plus propres?

        # Depot des pheromones
        for j in 1:length(C)
            if (cheminFourmis[idbestSol][j] ==1)
                #max(1,vecteurPheromones[j]* beta)
                vecteurPheromones[j] = min(1,vecteurPheromones[j]+ beta) # eviter que prob >1??

            end
            # Calcule les phéromones de chacune des fourmis
            #calculePheromones(vecteurPheromones, vecteurSol[j], bestSol)
            #vecteurPheromones[j] = MAJ_vectP(vecteurPheromones[j]) 
        end
        @show vecteurPheromones


        #@show vecteurDeToutesLesSolutionsFinales
        #@show vecteurPheromones

        #vecteurDeToutesLesSolutionsFinales = vecteurSolFourmis

        #vecteurSolFourmis[i] = vecteurPheromones[idbestSol]
        for j in eachindex(vecteurSolFourmis)
            if(vecteurSolFourmis[j] > meilleur)
                meilleur = vecteurSolFourmis[j] #meilleur ca doit etre int
            end
        end

        @show meilleur

    end

    # Recherche locale sur le meilleur trouvé au bout de n itérations 

    #resFinal = gloutonAmelioration(C, A, x, meilleur)
end



function roulette(nbAlea, vecteurSol)

    #val plafond(nbAlea/9)

    return ceil(nbAlea*length(vecteurSol))
    
end



function fourmisConstruction(vecteurSol, vecteurPheromones, C, A)

    vecteurSol::Vector{Int64}= zeros(length(vecteurPheromones))
    z::Int64 = 0
    lambda::Float64 = 5
    #i::Int64 = 1
    zTot::Int64 = 0

    vecteurSolPrime = Vector{String}()
    vecteurBool = Vector{Bool}(undef, length(vecteurPheromones))

    for i in 1:length(vecteurPheromones)
        push!(vecteurSolPrime,".")
    end

    taille::Int64 = length(vecteurPheromones)

    nbAlea = rand()

    posIndice1 = roulette(nbAlea, vecteurSol)
    #vecteurSol[Int(posIndice1)] = 1
    #println(vecteurSol)

    vecteurSol, z = fourmisConstruction2(C, A, vecteurSol, vecteurPheromones, Int(posIndice1))

    #println(vecteurSol)

    #=
    for i in 1:taille#eachindex(vecteurSol)    while (taille > 0) #
        println("\ni = $i")
        nbAlea = rand()
        
        if (vecteurSolPrime[i] == ".")
            
            if (nbAlea < vecteurPheromones[i])
                println("\ntruc")
                vecteurSol[i] = 1
                vecteurPheromones[i] = vecteurPheromones[i]*lambda
                println("\nAJOUT 1 -> $vecteurSol")


                println("\nvecteur solution -> $vecteurSol")
                vecteurSol, z = fourmisConstruction2(C, A, vecteurSol, vecteurPheromones)
                @show z

                zTot = zTot + z
                @show zTot
                vecteurPheromones = vecteurPheromones[setdiff(1:end, i)]

                place = findall(isequal(1), vecteurSol)
                vecteurSolPrime[place] .= "Occupé" 
                taille = taille -1
            else
                vecteurSol[i] = 0
                println("\nAJOUT 0 -> $vecteurSol")
                taille = taille -1
            end

        end

        #i = i+1

        println("\ntaille = $taille")
        
    end
=#
    return vecteurSol, z
    #=if (isAdmissible(C, A, vecteurSol) == true)
        return vecteurSol
    else
        return fourmisConstruction(vecteurSol, vecteurPheromones, C, A)
    end=#
end


function fourmisConstruction2(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, vecteurSol, vecteurPheromones, position)

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation
    n::Int64             = size(A,2)       # n nombre de variables
    m::Int64             = size(A,1)       # m nombre de contraintes

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    bestCandidate::Int64 = 0
    iteration::Int64 = 0


    # 1) REDUCTION DE L'INSTANCE SUR VARIABLES NON CONTRAINTES ----------------

    # Elimine de l'instance toutes les variables concernees par aucune contrainte
    variablesConcernees = findall(isequal(0), vec(sum(A, dims=1)))
    for j in variablesConcernees
        sol[j] = 1 # maj solution (partielle) en fixant a 1 la variable non contrainte
        z = z + C[j] # maj valeur de la solution (partielle)
    end

    # supprime les colonnes correspondant aux variables fixees
    index = index[setdiff(1:end, variablesConcernees)]
    C = C[setdiff(1:end, variablesConcernees)]
    A = A[:, setdiff(1:end, variablesConcernees)]

    # 2) CALCUL DES UTILITES --------------------------------------------------

   # candidates::Vector{Float64} = C ./ vec(sum(A, dims = 1))
    candidates = vecteurPheromones

    #println(" candidates -> $candidates")
    

    if verbose
        println("ivar   : ", index)
        #println("C      : ", C)
        #println("A      : ", A)
        println("U      : ", candidates)
        println(" ")
    end


    while (size(A,1) != 0) && (size(C,1) != 0) 
        iteration = iteration + 1

        # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Selection de l'indice du meilleur candidat au regard de son utilite
        if (iteration == 1)
            bestCandidate = position
        else
            bestCandidate = argmax(candidates)
        end
        #@show bestCandidate
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        #@show sol
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]
        #@show z
        # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

        # Identification des contraintes a supprimer du fait de la variable selectionnee
        lignestemp = findall(isequal(1), A[:,bestCandidate])

        #@show lignestemp

        # identifie toutes les colonnes qui doivent etre supprimées suite au candidat selectionne
        colonnetemp=(Int64)[]
        for i in lignestemp
            # scrute contrainte par contrainte les coefficients de A de valeur 1 (et elimine les eventuels doublons)
            colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
        end

        #@show colonnetemp
        # On supprime dans les structures index, A et C les valeurs correspondantes aux colonnes supprimées
        index      = index[setdiff(1:end, colonnetemp)]      # reduction du vecteur des indices des variables
        A          = A[:, setdiff(1:end, colonnetemp)]       # reduction des colonnes de la matrice des contraines
        C          = C[setdiff(1:end, colonnetemp)]          # reduction du vecteur de coefficients de la fonction objectif
        candidates = candidates[setdiff(1:end, colonnetemp)] # reduction du vecteur des utilites

        # On supprime dans la matrice A les lignes correspondantes aux contraintes impliquees par la variable selectionnee
        A          = A[setdiff(1:end, lignestemp), :]        # reduction des lignes de la matrice des contraines
        

        #@show A
        #@show C

        if verbose
            println("jselec : ", bestCandidate)
            println("ivar   : ", index)
            #println("C      : ", C)
            #println("A      : ", A)
            println("U      : ", candidates)
            println("-------- ")
        end
     
    end
    return sol, z
end


function isAdmissible(C, A, x)

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
        return false
        #@assert false "detection solution non-admissible"
    end
    #println( "admissible : oui | som(x_i) = ", length(var1), " ; z = ", z)
    println( "admissible")


    println( "\n $z")
    return true
end

#=
Si stagnation des resultats !COUP DE PIED!

function coupDePied()
end
=#