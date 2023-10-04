include("getfname.jl")

function glouton(C, A)

    #m lignes et n colonnes

    inutile, n = size(A) 

    (v1, v2, v3, v4, v5, v6, v7, v8, v9) = [view(A, :, i) for i in 1:size(A, 2)]

    V = [zeros(n) for _ in 1:n]
    V = (v1, v2, v3, v4, v5, v6, v7, v8, v9)

    coefDansLaFonctionObjective = update(A, C, V)

    supprimerLePlusGrosCoef(A, coefDansLaFonctionObjective)
#= 
    if (fonctionUtilite[1] == 0)
        println(" Tous les éléments ont été mis dans le sac")
    end
    else if {
        le sac est rempli sans avoir prit tous les éléments
    } 
    else {
        #mettre à jour la matrice A sinon boucle infinie
        glouton(C, A)
    }=#

end

#fonction qui retourne les occurences de chacunes des variables x1, x2, x3 ect..

function Occ(v, m, cmpt = 0)

    for i in 1:m
        cmpt = cmpt + v[i]
    end

    cmpt
end

# fonction qui retourne la fonction d'utilité en ordre décroissant

function calculFonctionUtilité(tableauOccurence, n)

    newtab::Vector{Float64} = zeros(Int, n)

    tOcc = merge_sort!(tableauOccurence)

    i::Int64 = 1
    j::Int64 = 9

    while i != 9
        newtab[i] = tOcc[j]
        j = j-1
        i = i+1
    end

    newtab
end


#TRI_FUSION de base

function merge_sort!(A, p = 1, r = length(A))
    if p < r
        q = div(p+r, 2)
        merge_sort!(A, p, q)
        merge_sort!(A, q+1, r)
        merge!(A, p, q, r)
    end
    A
end

function merge!(A, p, q, r)
    sentinel = typemax(eltype(A)) # une grand M
    L = A[p:q]
    R = A[(q+1):r]
    push!(L, sentinel)
    push!(R, sentinel)
    i, j = 1, 1
    for k in p:r
      if L[i] <= R[j]
          A[k] = L[i]
          i += 1
      else
          A[k] = R[j]
          j += 1
      end
    end
end

# Fonction qui retourne le coef le plus gros de la fonction objective 

function update(A, C, V)

    m, n = size(A)

    fonctionUtilite::Vector{Float64} = zeros(Int, n)
    tableauOccurence::Vector{Float64} = zeros(Int, n)
    cmpt::Int64 = 0


    for i in 1:n
        cmpt = Occ(V[i], m)
        tableauOccurence[i] = C[i]/cmpt
    end

    fonctionUtilite = calculFonctionUtilité(tableauOccurence, n)

    @show tableauOccurence
    @show fonctionUtilite

    println("\n")

    plusGrosCoef = fonctionUtilite[1]
    @show plusGrosCoef

    #recuperer l'indice du plus gros coef
    trouve::Bool = false
    i::Int64 = 0
    coefDansLaFonctionObjective::Int64 = 0

    while trouve != true

        i = i+1
        coefDansLaFonctionObjective = i

        if tableauOccurence[i] == plusGrosCoef
            trouve = true
        end

    end

    @show coefDansLaFonctionObjective
    #@show fonctionUtilite

    coefDansLaFonctionObjective
    
end


function supprimerLePlusGrosCoef(A, coefDansLaFonctionObjective)
    #=@show A
    m, n = size(A)
    APrime = zeros(Int, m, n)
    trouvé::Bool = false

    for i in 1:m
        for j in 1:n
            if(coefDansLaFonctionObjective == A[i,j])
                trouvé = true
            end
        end
    end



    @show trouvé

    for i in 1:n
        V[i][coefDansLaFonctionObjective]
    end=#

end