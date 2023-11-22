function ACO(C, A, nbIterationsACO, nbFourmis)

    # Initialisation --------------------------------------------

    vecteurSol::Vector{Int64}            = zeros(length(C))
    vecteurPheromones::Vector{Float64}   = Vector(undef, length(C)) #proba entre 0 et 1
    cheminFourmis::Vector{Vector{Int64}} = fill([], nbFourmis)      
    vecteurSolFourmis::Vector{Int64}     = zeros(nbFourmis)

    idbestSol::Int64   = 0 
    lambda::Float64    = 0.1  # coef evaporation des fourmis initié à 10%
    beta::Float64      = 1.1  # coef evaporation des fourmis
    meilleur::Float64  = 0
    z::Int64           = 0
    lastRestart::Int64 = 0

    #-------------------------------------------------------

    for i in 1:length(C)
        vecteurPheromones[i] = 1/length(C)
    end

    #-----------------------------------------------------------

    # Boucle POUR pour avoir n itérations ------------------------
    for i in 1:nbIterationsACO 

        idbestSol = 0

        # Boucle POUR utilisée pour chaque fourmis --------------
        for j in 1:nbFourmis
            cheminFourmis[j], z = fourmisConstruction(vecteurSol, vecteurPheromones, C, A)
            vecteurSolFourmis[j] = sum(cheminFourmis[j][m]*C[m] for m in 1:length(C))
        end

        #-------------------------------------------------------

        idbestSol = argmax(vecteurSolFourmis) 

        #-------------------------------------------------------

        # Recherche locale sur le meilleur trouvé au bout de n itérations 

        cheminFourmis[idbestSol], inutile = gloutonAmelioration(C, A, cheminFourmis[idbestSol], meilleur)


        #-------------------------------------------------------

        vecteurPheromones = coupDePied(vecteurPheromones, cheminFourmis, idbestSol, i, nbIterationsACO, lastRestart, lambda, beta, meilleur)
       
        #-------------------------------------------------------

        for j in eachindex(vecteurSolFourmis)
            if(vecteurSolFourmis[j] > meilleur)
                meilleur = vecteurSolFourmis[j]
            end
        end

        @show meilleur

    end
end


#Si stagnation des resultats !COUP DE PIED!

function coupDePied(vecteurPheromones, cheminFourmis, idbestSol, iter, iterMax, lastRestart, lambda, beta, meilleur)

    # Evaporation ------------------------------------------
    for j in 1:length(vecteurPheromones)
        vecteurPheromones[j] = vecteurPheromones[j] * lambda    
    end

    #-------------------------------------------------------

    # Depot des pheromones ---------------------------------
    for j in 1:length(vecteurPheromones)
        if (cheminFourmis[idbestSol][j] ==1)
            
            vecteurPheromones[j] = min(1,vecteurPheromones[j]+ beta) # eviter que prob >1??
        end
    end

    #-------------------------------------------------------

    # Coup de pied -----------------------------------------
    #si solution stagnante #si pheromones a 0
    #si on peut donner un coup de pied (comparer nb iteration avec nb iteration max et le dernier coup de pieds)

    if (meilleur == cheminFourmis[idbestSol]) && (findall(isequal(0), vecteurPheromones) == true) && (iter-lastRestart > 10)
        println(" COUP DE PIED!!!")

        # Perturbation 1 -----------------------------------
        for j in 1:length(vecteurPheromones)
            vecteurPheromones[j] = vecteurPheromones[j]*0.95*(log10(iter)/log10(iterMax)) # POUR TOUT CASSERRRR
        end
        #-------------------------------------------------------

        # Perturbation 2 ------------------------------------

        for j in rand(0:length(vecteurPheromones))
            vecteurPheromones[rand(1:length(vecteurPheromones))] = rand(.05:.1:(iter-(1/iterMax))*.5, 1, 1)
         end

        for j in 1:length(vecteurPheromones)
            if (vecteurPheromones[j] < 0.1) 
                vecteurPheromones[j] =  rand(.05:.1:(iter-(1/iterMax))*.5, 1, 1)
            end
        end

        #-------------------------------------------------------

    end
    return vecteurPheromones
end



function roulette(nbAlea, vecteurSol)
    return ceil(nbAlea*length(vecteurSol))
end



function fourmisConstruction(vecteurSol, vecteurPheromones, C, A)

    z::Int64        = 0
    nbAlea::Float64 = rand()

    #-------------------------------------------------------
    
    posIndice1    = roulette(nbAlea, vecteurSol)
    vecteurSol, z = fourmisConstructionGlouton(C, A, vecteurSol, vecteurPheromones, Int(posIndice1))

    return vecteurSol, z
end


function fourmisConstructionGlouton(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, vecteurSol, vecteurPheromones, position)

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
