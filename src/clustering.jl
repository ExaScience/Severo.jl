# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, nonzeros, getcolptr, rowvals, nnz
import StatsBase

include("modularity.jl")

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
@partial cluster(SNN::NeighbourGraph; kw...) = cluster(default_rng(), SNN; kw...)

"""
    cluster(rng::AbstractRNG, SNN::NeighbourGraph; algorithm=:louvain, resolution=0.8, nstarts=1, niterations=10) where T

Cluster cells based on a neighbourhood graph.

**Arguments**:

    - `rng`: random number generator
    - `SNN`: shared neighbours graph
    - `algorithm`: clustering algorithm to use (louvain)
    - `resolution`: parameters above 1 will lead to larger communities whereas below 1 lead to smaller ones
    - `nstarts`: number of random starts
    - `niterations`: maximum number of iterations per random start
    - `group_singletons`: group singletons into nearest cluster, if false keeps singletons

**Return values**:

cluster assignment per cell
"""
function cluster(rng::AbstractRNG, SNN::NeighbourGraph; algorithm=:louvain, resolution=0.8, nrandomstarts=10, niterations=10, verbose=false, group_singletons::Bool=true)
    if isa(algorithm, AbstractString)
        algorithm = Symbol(algorithm)
    end

    assignment = if algorithm == :louvain
        n = Network(SNN.array)
        cl = louvain_multi(rng, n, nrandomstarts, resolution=resolution, max_iterations=niterations, verbose=verbose)
        cl.nodecluster
    else
        error("unknown clustering algorithm: $algorithm")
    end

    if group_singletons
        group_singletons!(assignment, SNN.array)
    end

    NamedArray(assignment, (SNN.dicts[1],), (:cells,))
end

function group_singletons!(lbls::Vector{T}, snn::SparseMatrixCSC{R, Int64}; min_count::Int=1) where {T <: Integer, R <: Real}
    @assert length(lbls) == size(snn,1)

    c = StatsBase.counts(lbls)
    singletons = Set(findall(c .<= min_count))

    if !isempty(singletons)
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

        relabel!(lbls, length(c))
    end

    lbls
end

function group_singletons!(lbls::NamedVector{T}, SNN::NeighbourGraph; min_count::Int=1) where {T <: Integer}
    group_singletons!(lbls.array, SNN.array)
    lbls
end

group_singletons(lbls::NamedVector{T}, SNN::NeighbourGraph; min_count::Int=1) where {T <: Integer} = group_singletons!(copy(lbls), SNN; min_count=min_count)
