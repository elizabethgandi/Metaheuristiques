function destruction(C, xEntree, zEntree)

    x = copy(xEntree)
    z = zEntree

    Liste1 = findall(x->x==1, x)         # Identifie les indices de tous les x[i]=1
    nListe1 = length(Liste1)             # Nombre de 1 dans la solution

    shuffle!(Liste1)                     # Melange la liste d'indices
    nDetruits = rand(1:nListe1)          # Tire aleatoirement le nombre de 1 a detruire
    ListeDetruits = Liste1[1:nDetruits]  # Extrait la liste d'indices a detruire

    # MISE A JOUR DE LA SOLUTION
    for i in ListeDetruits
        x[i] = 0
        z -= C[i]
    end
    return x, z
end
