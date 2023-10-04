include("getfname.jl")

function glouton(C, A)

    #m lignes et n colonnes

    m, n = size(A) 

    x::Vector{Int64}= zeros(Int, n)
    fonctionUtilite::Vector{Float64} = zeros(Int, n)
    tableauOccurence::Vector{Float64} = zeros(Int, n)

    (v1, v2, v3, v4, v5, v6, v7, v8, v9) = [view(A, :, i) for i in 1:size(A, 2)]

    V = [zeros(n) for _ in 1:n]
    V = (v1, v2, v3, v4, v5, v6, v7, v8, v9)

    for i in 1:n
        cmpt = Occ(V[i], m)
        tableauOccurence[i] = C[i]/cmpt
    end

    #@show tableauOccurence

    fonctionUtilite = calculFonctionUtilité(tableauOccurence, n)

    @show fonctionUtilite

end

function Occ(v, m)
    cmpt = 0

    for i in 1:m
        cmpt = cmpt + v[i]
    end

    return cmpt
end


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


#TRI_FUSION

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
    end
end


#= croissant
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
    end
end=#