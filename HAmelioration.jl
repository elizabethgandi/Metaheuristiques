function gloutonAmelioration(C, A, xConstruction, zConstruction)
    tabDesIndicesVariablesAZero = findall(isequal(0), xConstruction)
    tabDesIndicesVariablesAUn   = findall(isequal(1), xConstruction)

    verbose  = true
    zCourant = 0
    contraintesSaturees = zeros(Inf64, size(A,1))
    zOptimal = zConstruction

    # 2-1 echange --------------------------------------------------

    if (verbose)
        println("> 2-1 : ")
    end
    for i in tabDesIndicesVariablesAUn
        for j in tabDesIndicesVariablesAUn 
            if (i<j)
                for k in tabDesIndicesVariablesAZero

                    # Tester si la solution zCourant est meilleure que zBest
                    if (zCourant-C[i]-C[j]+C[k] > zOptimal)
                        #Remplissage a la main des contraintes saturees ET verification d'admissibilite
                        if ((findfirst(x->x>1,contraintesSaturees - A[:,i] - A[:,j] + A[:,k])) == nothing) #&& isAdmissible(xAmelioration) == true)
                            #v10A
                            v1 = i
                            #v10b
                            v2 = j
                            #v01
                            v3 = k
                        end
                    end
                end 
            end
        end
    end


end

#function isAdmissible(xSolution)
#
#end

function isAd(C, A, x)

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
        @assert false "detection solution non-admissible"
    end
    println( "admissible : oui | som(x_i) = ", length(var1), " ; z = ", z)
    return true

end

function isAllowed(A::SparseMatrixCSC{Int64, Int64}, x::SparseVector{Int64})
    index = findall(!iszero, x)
    constraints = Vector{Int64}(undef, 0)
    # Pour chaque variable, regarder dans quelles contraintes elle apparait
    for i in index
        append!(constraints, findall(t->t==1, A[:,i]))
    end
    unique!(constraints)
    # Verifie chacune des contraintes
    for i in constraints
        somme::Int64 = 0
        for j in index
            somme = somme + A[i,j]
        end
        if somme > 1
            return false
        end
    end
    return true
end