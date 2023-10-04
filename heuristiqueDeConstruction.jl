include("getfname.jl")

function glouton(C, A)

    inutile, n = size(A) 
    @show n

    x::Vector{Int64}= zeros(Int, n)
    fonctionUtilite::Vector{Int64} = zeros(Int, n)
    tableauOccurence::Vector{Int64} = zeros(Int, n)

    @show x
    @show fonctionUtilite
    @show tableauOccurence


    #=@show x
    @show fonctionUtilite
    @show tableauOccurence=#


end