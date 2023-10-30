# --------------------------------------------------------------------------- #
# collect the un-hidden filenames available in a given directory

function getfname(target)
  # target : string := chemin + nom du repertoire ou se trouve les instances
  pathCourant = pwd()

  # positionne le currentdirectory dans le repertoire cible
  cd(target)

  # retourne le repertoire courant
  #println("pwd = ", pwd())

  # recupere tous les fichiers se trouvant dans le repertoire data
  allfiles = readdir()

  # vecteur booleen qui marque les noms de fichiers valides
  flag = trues(size(allfiles))

  k=1
  for f in allfiles
      # traite chaque fichier du repertoire
      if f[1] != '.'
          # pas un fichier cache => conserver
          println("fname = ", f)
      else
          # fichier cache => supprimer
          flag[k] = false
      end
      k = k+1
  end

  cd(pathCourant)

  # extrait les noms valides et retourne le vecteur correspondant
  finstances = allfiles[flag]
  return finstances
end

# --------------------------------------------------------------------------- #
# Loading an instance of SPP (format: OR-library)

function loadSPP(fname)
    f=open(fname)
    # lecture du nbre de contraintes (m) et de variables (n)
    m, n = parse.(Int, split(readline(f)) )
    # lecture des n coefficients de la fonction economique et cree le vecteur d'entiers c
    C = parse.(Int, split(readline(f)) )
    # lecture des m contraintes et reconstruction de la matrice binaire A
    A=zeros(Int, m, n)
    for i=1:m
        # lecture du nombre d'elements non nuls sur la contrainte i (non utilise)
        readline(f)
        # lecture des indices des elements non nuls sur la contrainte i
        for valeur in split(readline(f))
          j = parse(Int, valeur)
          A[i,j]=1
        end
    end
    close(f)
    return C, A
end
