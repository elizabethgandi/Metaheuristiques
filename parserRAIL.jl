# ---- Parser reading an instance "rail" from OR-library 
#      (http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)

function loadInstanceRAIL(fname::String)

    f = open(fname)
    nbctr, nbvar = parse.(Int, split(readline(f)))
    A = zeros(Int, nbctr, nbvar)       # matrix of constraints
    c = zeros(Int, nbvar)              # vector of costs
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                A[j,i] = 1
            end
        end
    end
    close(f)
    return nbvar, nbctr, A, c
end