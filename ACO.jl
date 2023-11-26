function ACO(C, A, nbIterationsACO, nbFourmis)

    # Initialisation ----------------------------------------------------------

    vecteurSol::Vector{Int64}            = zeros(length(C))          # Vecteur de solution pour une fourmis
    vecteurPheromones::Vector{Float64}   = Vector(undef, length(C))  # Vecteur de phéromones
    cheminFourmis::Vector{Vector{Int64}} = fill([], nbFourmis)       # Chemin emprunté par une fourmis
    vecteurSolFourmis::Vector{Int64}     = zeros(nbFourmis)          # Vecteur contenant les solutions de toutes les foumis pour 1 itération

    idbestSol::Int64   = 0    # Indice de la meilleur solution trouvée par une fourmis à une itération i
    lambda::Float64    = 0.2  # Coefficient d' évaporation des fourmis initié à 10%
    beta::Float64      = 0.7  # Coefficient de dépots des fourmis sur le vecteur phéromone
    meilleur::Float64  = 0    # Meilleur solution z pour l'itération i
    z::Int64           = 0    # Valeur de la fonction objective
    lastRestart::Int64 = 0    # Valeur pour le coup de pied

    verbose::Bool      = false

    #--------------------------------------------------------------------------
    # Plot --------------------------------------------------------------------

    x         = [] # Nombres d' itérations en ordonnées
    y         = [] # La valeur de 1 fourmis en abscisse
    yMeilleur = [] # La meilleure des fourmis à 1 itération en abscisse
    xMeilleur = [] # L'itération correspondant à la meilleur des solutions trouvées

    #--------------------------------------------------------------------------
    # Initialisation du vecteur phéromone -------------------------------------

    for i in 1:length(C)
        vecteurPheromones[i] = 1/length(C)
    end

    #--------------------------------------------------------------------------
    # Boucle correspondant aux i itérations -----------------------------------
    
    for i in 1:nbIterationsACO 

        idbestSol = 0

        print("\niter $i -> z = $z \n")

        #----------------------------------------------------------------------
        # Boucle correspondant au lancement de chacunes des fourmis -----------

        for j in 1:nbFourmis

            # Construction du chemin emprunté par les fourmis
            cheminFourmis[j], z  = fourmisConstruction(C, A, vecteurSol, vecteurPheromones)

            # Stockage de la valeur de la fonction objective z dans le vecteur de solutions des fourmis
            vecteurSolFourmis[j] = z 

            #------------------------------------------------------------------
            # Insertion des valeurs de z et i pour l'affichage de plot --------
            push!(y,z)
            push!(x,i)

        end

        #----------------------------------------------------------------------
        # Identification de la meilleur solution trouvée parmis toutes les fourmis
        idbestSol = argmax(vecteurSolFourmis) 

        #----------------------------------------------------------------------
        # Gestion du vecteur de phéromones ------------------------------------
        vecteurPheromones = gestionEtCoupDePied(vecteurPheromones, cheminFourmis, idbestSol, i, nbIterationsACO, lastRestart, lambda, beta, meilleur, z)
       
        #----------------------------------------------------------------------
        # Identification du meilleur z trouvé depuis le lancement du programme

        for j in eachindex(vecteurSolFourmis)
            if (vecteurSolFourmis[j] > meilleur)
                meilleur = vecteurSolFourmis[j]
            end
        end

        #------------------------------------------------------------------
        # Insertion des valeurs de zMeilleur et i pour l'affichage de plot 

        push!(yMeilleur, meilleur)
        push!(xMeilleur, i)

        if verbose
            #println(" cheminFourmis   : ", cheminFourmis)
            println(" zBest      : ", meilleur)
            println(" lancé      : ", i)
            println(" ")
        end

    end

    #--------------------------------------------------------------------------
    # Affichage plot ----------------------------------------------------------

    scatter(x,y, s=2, c="black")
    plt.plot(xMeilleur,yMeilleur, c="red")
    plt.show()

    #--------------------------------------------------------------------------
    # Activation de notre recherche locale sur le meilleur des résultats trouvés

    cheminFourmis[idbestSol], inutile = gloutonAmelioration(C, A, cheminFourmis[idbestSol], meilleur)

    return cheminFourmis[idbestSol], meilleur
end


function gestionEtCoupDePied(vecteurPheromones, cheminFourmis, idbestSol, iter, iterMax, lastRestart, lambda, beta, meilleur, zFourmis)

    #--------------------------------------------------------------------------
    # Evaporation -------------------------------------------------------------

    for j in 1:length(vecteurPheromones)
        vecteurPheromones[j] = vecteurPheromones[j] * lambda    
    end

    #--------------------------------------------------------------------------
    # Dépôt des pheromones ----------------------------------------------------

    for j in 1:length(vecteurPheromones)
        if (cheminFourmis[idbestSol][j] ==1)
            vecteurPheromones[j] = min(1,vecteurPheromones[j]+ beta) 
        end
    end

    #--------------------------------------------------------------------------
    # Coup de pied si stagnation ----------------------------------------------

    if (meilleur == zFourmis) && ((length((findall(isequal(1.0), vecteurPheromones)))>0) == true) && ((mod(iter-lastRestart, 10))==0)
        #println(" COUP DE PIED!!!")

        #----------------------------------------------------------------------
        # Perturbation 1 ------------------------------------------------------
        for j in 1:length(vecteurPheromones)
            vecteurPheromones[j] = vecteurPheromones[j]*0.95*(log10(iter)/log10(iterMax)) # POUR TOUT CASSERRRR
        end

        #----------------------------------------------------------------------
        # Perturbation 2 ------------------------------------------------------

        for j in rand(0:length(vecteurPheromones))
            vecteurPheromones[rand(1:length(vecteurPheromones))] = rand(.05:(iter-(1/iterMax))*.5)
        end

        for j in 1:length(vecteurPheromones)
            if (vecteurPheromones[j] < 0.1) 
                vecteurPheromones[j] =  rand(.05:(iter-(1/iterMax))*.5)
            end
        end

        #----------------------------------------------------------------------

    end
    return vecteurPheromones
end



function roulette(nbAlea, vecteurSol)

    return ceil(nbAlea*length(vecteurSol))

end


function positionRoulette(vecteurPheromones)

    cumule::Vector{Float64} = cumsum(vecteurPheromones) 
    vTiree::Float64         = cumule[end]*rand()
    i::Int64                = 1

    while (cumule[i] < vTiree)
        i = i+1
    end

    return i
end


function fourmisConstruction(C_entree::Vector{Int64}, A_entree::Matrix{Int64}, vecteurSol, vecteurPheromones)

    C = copy(C_entree) 
    A = copy(A_entree) 

    verbose::Bool = false

    # Initialisation ----------------------------------------------------------

    n::Int64             = size(A,2)       # n nombre de variables

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    bestCandidate::Int64 = 0
    iteration::Int64     = 0

    # -------------------------------------------------------------------------
    # ACO ---------------------------------------------------------------------

    nbAlea::Float64   = rand()
    position::Float64 = roulette(nbAlea, vecteurSol)

    # -------------------------------------------------------------------------

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
    candidates = copy(vecteurPheromones)

    while (size(A,1) != 0) && (size(C,1) != 0) 
        iteration = iteration + 1

        # 3) CHOIX DU CANDIDAT A AJOUTER A LA SOLUTION COURANTE ---------------

        # Sélection de l'indice du candidat, soit à une certaine position "postion" donnée
        # pour la première itération ou dans la roulette sinon
        if (iteration == 1)
            bestCandidate = Int(position)
        else
            #bestCandidate = argmax(candidates)
            bestCandidate = positionRoulette(candidates)
        end
        # Mise à jour de la solution avec le candidat selectionne
        sol[index[bestCandidate]] = 1
        # Mise à jour de la valeur de la fonction objective avec le candidat selectionne
        z = z + C[bestCandidate]
        # 4) REDUCTION DU PROBLEME SUITE AU CANDIDAT SELECTIONNE --------------

        # Identification des contraintes a supprimer du fait de la variable selectionnee
        lignestemp = findall(isequal(1), A[:,bestCandidate])

        # identifie toutes les colonnes qui doivent etre supprimées suite au candidat selectionne
        colonnetemp=(Int64)[]
        for i in lignestemp
            # scrute contrainte par contrainte les coefficients de A de valeur 1 (et elimine les eventuels doublons)
            colonnetemp = union(colonnetemp, findall(isequal(1), A[i,:]))
        end

        # On supprime dans les structures index, A et C les valeurs correspondantes aux colonnes supprimées
        index      = index[setdiff(1:end, colonnetemp)]      # reduction du vecteur des indices des variables
        A          = A[:, setdiff(1:end, colonnetemp)]       # reduction des colonnes de la matrice des contraines
        C          = C[setdiff(1:end, colonnetemp)]          # reduction du vecteur de coefficients de la fonction objectif
        candidates = candidates[setdiff(1:end, colonnetemp)] # reduction du vecteur des utilites

        # On supprime dans la matrice A les lignes correspondantes aux contraintes impliquees par la variable selectionnee
        A          = A[setdiff(1:end, lignestemp), :]        # reduction des lignes de la matrice des contraines
        
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
