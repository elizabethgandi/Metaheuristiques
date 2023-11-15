function ACO(C, A, nbIterationsACO, nbFourmis)
    # Initialisation

    vecteurSol::Vector{Int64}          = zeros(length(C))
    vecteurPropa::Vector{Int64}        = zeros(length(C))
    vecteurPheromones::Vector{Float64} = Vector(undef, length(C))
    idbestSol::Int64                   = 0 
    lambda::Float64                    = 0.9 # coef evaporation des fourmis
    beta::Float64                      = 1.1  # coef evaporation des fourmis
    meilleur::Float64 = 0

    cheminFourmis::Vector{Vector{Int64}} = fill([], length(C))
    vecteurSolFourmis::Vector{Int64}     = zeros(length(C))
    vecteurDeToutesLesSolutionsFinales::Vector{Float64} = zeros(length(C))

    
    print("ca passe dedans")

    for i in 1:length(C)
        vecteurPheromones[i] = 1/2
    end

    # Premier appel pour obtenir une fonction d'utilité des phéromones
    #vecteurSol = gloutonConstruction(C, A)

    # Boucle POUR pour avoir n itérations
    for i in 1:nbIterationsACO # fixer nbIterationsACO à 1 pour l'instant?
        idbestSol = 0
        # Boucle POUR utilisée pour chaque fourmis
        for j in 1:nbFourmis

            # Fonction qui construit une solution avec la roulette
            cheminFourmis[j] = fourmisConstruction(vecteurSol, vecteurPheromones)
            
            vecteurSolFourmis[j] = sum(cheminFourmis[j][m]*C[m] for m in 1:length(C)) #a test la syntaxe
       
        end

        idbestSol = argmax(vecteurSolFourmis) 

       

        # Evaporation
        for j in 1:length(vecteurPheromones)
            vecteurPheromones[j] = vecteurPheromones[j] * lambda
        end


        # Depot des pheromones
        for j in 1:length(cheminFourmis)
            if (cheminFourmis[j] ==1)
                max(1,vecteurPheromones[j]* beta)
                vecteurPheromones[j] = vecteurPheromones[j]* beta

            end
            # Calcule les phéromones de chacune des fourmis
            #calculePheromones(vecteurPheromones, vecteurSol[j], bestSol)
            #vecteurPheromones[j] = MAJ_vectP(vecteurPheromones[j]) 
        end


        @show vecteurDeToutesLesSolutionsFinales
        @show vecteurPheromones


        vecteurDeToutesLesSolutionsFinales[i] = vecteurPheromones[idbestSol]

        if(vecteurDeToutesLesSolutionsFinales[i] > meilleur)
            meilleur = vecteurDeToutesLesSolutionsFinales[i] #meilleur ca doit etre int
        end

    end
end

function fourmisConstruction(vecteurSol, vecteurPheromones)

    # ajouter le brut force des meilleures vecteurDeToutesLesSolutionsFinales

    for i in eachindex(vecteurSol)

        nbAlea = rand()
        
        if (nbAlea < vecteurPheromones[i])
            vecteurSol[i] = 1
        else
            vecteurSol[i] = 0
        end
    end
    
    return vecteurSol
end

#=
Si stagnation des resultats !COUP DE PIED!

function coupDePied()
end
=#