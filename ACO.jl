function ACO(C, A, nbIterationsACO, nbFourmis)

    # Initialisation --------------------------------------------

    vecteurSol::Vector{Int64}            = zeros(length(C))
    vecteurPheromones::Vector{Float64}   = Vector(undef, length(C)) #proba entre 0 et 1
    cheminFourmis::Vector{Vector{Int64}} = fill([], nbFourmis)      
    vecteurSolFourmis::Vector{Int64}     = zeros(nbFourmis)

    idbestSol::Int64   = 0 
    lambda::Float64    = 0.2  # coef evaporation des fourmis initié à 10%
    beta::Float64      = 0.7  # coef de dépots des fourmis
    meilleur::Float64  = 0
    z::Int64           = 0
    lastRestart::Int64 = 0
    verbose::Bool      = false

    #-------------------------------------------------------
    # Plot -------------------------------------------------

    x         = [] # itérations
    y         = [] # la valeur de 1 fourmis
    yMeilleur = [] # la meilleure des fourmis à 1 itération
    xMeilleur = []

    #-------------------------------------------------------

    for i in 1:length(C)
        vecteurPheromones[i] = 1/length(C)
    end

    #-----------------------------------------------------------

    # Boucle POUR pour avoir n itérations ------------------------
    for i in 1:nbIterationsACO 

        idbestSol = 0

        print("\niter $i -> z = $z")


        # Boucle POUR utilisée pour chaque fourmis --------------
        for j in 1:nbFourmis
            cheminFourmis[j], z = fourmisConstruction(vecteurSol, vecteurPheromones, C, A)
            push!(y,z)
            push!(x,i)
            vecteurSolFourmis[j] = sum(cheminFourmis[j][m]*C[m] for m in 1:length(C)) # remplacer par z?
        end

        #-------------------------------------------------------

        idbestSol = argmax(vecteurSolFourmis) 

        #-------------------------------------------------------

        # Recherche locale sur le meilleur trouvé au bout de n itérations 

        # Amélioration trop longue!!!!!

        #cheminFourmis[idbestSol], inutile = gloutonAmelioration(C, A, cheminFourmis[idbestSol], meilleur)
        @show cheminFourmis[idbestSol]


        #-------------------------------------------------------

        vecteurPheromones = coupDePied(vecteurPheromones, cheminFourmis, idbestSol, i, nbIterationsACO, lastRestart, lambda, beta, meilleur, z)
       
        #-------------------------------------------------------


        for j in eachindex(vecteurSolFourmis)
            if (vecteurSolFourmis[j] > meilleur)
                meilleur = vecteurSolFourmis[j]
            end
        end

        push!(yMeilleur, meilleur)
        push!(xMeilleur, i)
        #@show yMeilleur

        if verbose
            #println(" cheminFourmis   : ", cheminFourmis)
            println(" zBest      : ", meilleur)
            println(" lancé      : ", i)
            println(" ")
        end

       

    end

    #x = collect(1:length(y))
    #@show y
    #plot(x,y)
    scatter(x,y, s=2, c="black")
    scatter(xMeilleur,yMeilleur, s=4, c="red")

    #cheminFourmis[idbestSol], inutile = gloutonAmelioration(C, A, cheminFourmis[idbestSol], meilleur)

    return cheminFourmis[idbestSol], meilleur
end


function coupDePied(vecteurPheromones, cheminFourmis, idbestSol, iter, iterMax, lastRestart, lambda, beta, meilleur, zFourmis)

    # Evaporation ------------------------------------------
    for j in 1:length(vecteurPheromones)
        vecteurPheromones[j] = vecteurPheromones[j] * lambda    
    end

    #-------------------------------------------------------

    # Dépôt des pheromones ---------------------------------
    for j in 1:length(vecteurPheromones)
        if (cheminFourmis[idbestSol][j] ==1)
            vecteurPheromones[j] = min(1,vecteurPheromones[j]+ beta) #+
        end
    end

    #-------------------------------------------------------

    if (meilleur == zFourmis) && ((length((findall(isequal(1.0), vecteurPheromones)))>0) == true) && ((mod(iter-lastRestart, 10))==0)
        #println(" COUP DE PIED!!!")

        # Perturbation 1 -----------------------------------
        for j in 1:length(vecteurPheromones)
            vecteurPheromones[j] = vecteurPheromones[j]*0.95*(log10(iter)/log10(iterMax)) # POUR TOUT CASSERRRR
        end
        #-------------------------------------------------------

        # Perturbation 2 ------------------------------------

        for j in rand(0:length(vecteurPheromones))
            vecteurPheromones[rand(1:length(vecteurPheromones))] = rand(.05:(iter-(1/iterMax))*.5)
        end

        for j in 1:length(vecteurPheromones)
            if (vecteurPheromones[j] < 0.1) 
                vecteurPheromones[j] =  rand(.05:(iter-(1/iterMax))*.5)
            end
        end

        #-------------------------------------------------------

    end
    return vecteurPheromones
end



function roulette(nbAlea, vecteurSol)
    return ceil(nbAlea*length(vecteurSol))
end

function positionRoulette(vecteurPheromones)
    cumule = cumsum(vecteurPheromones) 
    vTiree = cumule[end]*rand()
    i = 1
    while (cumule[i] < vTiree)
        i = i+1
    end
    return i
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

    index::Vector{Int64} = collect(1:n)    # Index d'origines des variables
    sol::Vector{Int64}   = zeros(Int64,n)  # Vecteur de base de la solution
    z::Int64             = 0               # z la valeur de la fonction objective

    bestCandidate::Int64 = 0
    iteration::Int64     = 0


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

        # Selection de l'indice du meilleur candidat au regard de son utilite
        if (iteration == 1)
            bestCandidate = position
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
