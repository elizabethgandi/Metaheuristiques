include("getfname.jl")

function glouton(C, A)

    #m ligne et n colone

    m, n = size(A) 
    @show n

    x::Vector{Int64}= zeros(Int, n)
    fonctionUtilite::Vector{Int64} = zeros(Int, n)
    tableauOccurence::Vector{Float64} = zeros(Int, n)

    @show x
    @show fonctionUtilite
    @show tableauOccurence


    #=@show x
    @show fonctionUtilite
    @show tableauOccurence=#
    (v1, v2, v3, v4, v5, v6, v7, v8, v9) = [view(A, :, i) for i in 1:size(A, 2)]

    @show v1
    @show v2
    @show v3
    @show v4
    @show v5
    @show v6
    @show v7
    @show v8
    @show v9


    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v1[i]
    end

    tableauOccurence[1] = C[1]/cmpt

    #

    @show cmpt

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v2[i]
    end

    tableauOccurence[2] = C[2]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v3[i]
    end

    tableauOccurence[3] = C[3]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v4[i]
    end

    tableauOccurence[4] = C[4]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v5[i]
    end

    tableauOccurence[5] = C[5]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v6[i]
    end

    tableauOccurence[6] = C[6]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v7[i]
    end

    tableauOccurence[7] = C[7]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v8[i]
    end

    tableauOccurence[8] = C[8]/cmpt

    @show cmpt

    #

    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v9[i]
    end

    tableauOccurence[9] = C[9]/cmpt

    @show cmpt

    #
    @show tableauOccurence

    print("\n")
    
    print("\n")
    
    print("\n")

    cp1 = Occ(v1, m)
    cp2 = Occ(v9, m)
    
end


function tabOccurences(tableauOccurence, m, cmpt)
    #cmpt = 0

    for i in 1:m
        tableauOccurence[i] = C[i]/cmpt
    end

    @show cmpt
    @show tableauOccurence

    return tableauOccurence
end

function Occ(v, m)
    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v[i]
    end

    @show cmpt

    return cmpt
end

#=

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
    sentinel = typemax(eltype(A))
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
    end=#
#end