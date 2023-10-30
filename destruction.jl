function destruction(C, xEntree, zEntree)

    x = copy(xEntree)
    z = zEntree

    Liste1 = findall(x->x==1, x)  # identifie les indices de tous les x[i]=1
    nListe1 = length(Liste1)      # nombre de 1 dans la solution

    shuffle!(Liste1)              # melange la liste d'indices
    nDetruits = rand(1:nListe1)   # tire aleatoirement le nombre de 1 a detruire
    ListeDetruits = Liste1[1:nDetruits] # extrait la liste d'indices a detruire
    #println("detruit $nDetruits variables à 1 en position(s) ",ListeDetruits)

    # met a jour la solution
    for i in ListeDetruits
        x[i] = 0
        z -= C[i]
    end
    return x, z
end
