# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, nonzeros, getcolptr, rowvals, nnz
import StatsBase

include("modularity.jl")

function modularity_cluster(SNN::SparseMatrixCSC{Float64, Int64}; modularity=1, resolution=0.8, algorithm=1,
        randomseed=0, nrandomstarts=10, niterations=10, verbose=true)
    m,n = size(SNN)

    assignment = Vector{Int32}(undef, max(m,n));
    ccall(("ModularityClustering", libcell), Cvoid,
        (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Int32, Int32, Int64,
        Int32, Float64, Int32, Int32, Int32, Int32, Bool, Ptr{Int32}),
        getcolptr(SNN), rowvals(SNN), nonzeros(SNN), m, n, nnz(SNN),
        modularity, resolution, algorithm, nrandomstarts, niterations,
        randomseed, verbose, assignment)

    Int64.(assignment .+ 1)
end

"""
    cluster(SNN::NeighbourGraph; algorithm=:louvain, resolution=0.8, nstarts=1, niterations=10) where T

Cluster cells based on a neighbourhood graph.

**Arguments**:

    - `SNN`: shared neighbours graph
    - `algorithm`: clustering algorithm to use (louvain)
    - `resolution`: parameters above 1 will lead to larger communities whereas below 1 lead to smaller ones
    - `nstarts`: number of random starts
    - `niterations`: maximum number of iterations per random start
    - `group_singletons`: group singletons into nearest cluster, if false keeps singletons

**Return values**:

cluster assignment per cell
"""
function cluster(SNN::NeighbourGraph; algorithm=:louvain, resolution=0.8, nrandomstarts=10, niterations=10, verbose=false, group_singletons::Bool=true)
    if isa(algorithm, AbstractString)
        algorithm = Symbol(algorithm)
    end

    assignment = if algorithm == :louvaincpp
        modularity_cluster(SNN.array, algorithm=1, resolution=resolution,
                nrandomstarts=nrandomstarts, niterations=niterations, verbose=verbose)
    elseif algorithm == :louvain
        n = Network(SNN.array)
        cl = louvain_multi(n, nrandomstarts, resolution=resolution, max_iterations=niterations, verbose=verbose)
        cl.nodecluster
    else
        error("unknown clustering algorithm: $algorithm")
    end

    if group_singletons
        group_singletons!(assignment, SNN)
    end

    NamedArray(assignment, (SNN.dicts[1],), (:cells,))
end

function group_singletons!(lbls::NamedVector{T}, SNN::NeighbourGraph; min_count::Int=1) where {T <: Integer}
    @assert length(lbls) == size(SNN,1)
    snn = SNN.array

    c = StatsBase.counts(lbls)
    singletons = Set(findall(c .<= min_count))
    others = findall(c .> min_count)

    singles = findall(in(singletons), lbls)
    for j in singles
        val, closest = findmax(others) do cl
            nzv = nonzeros(snn)
            rv = rowvals(snn)
            neighbours = findall(idx -> lbls[rv[idx]] == cl, nzrange(snn, j))
            mean(nzv[neighbours]) # XXX does not make sense to me
        end

        lbls[j] = closest
    end

    lbls
end

group_singletons(lbls::NamedVector{T}, SNN::NeighbourGraph; min_count::Int=1) where {T <: Integer} = group_singletons!(copy(lbls), SNN; min_count=min_count)
