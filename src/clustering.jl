import SparseArrays: SparseMatrixCSC, nonzeros, getcolptr, rowvals, nnz

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

    assignment .+= 1
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

**Return values**:

cluster assignment per cell
"""
function cluster(SNN::NeighbourGraph; algorithm=:louvain, resolution=0.8, nrandomstarts=1, niterations=10)
    if isa(algorithm, AbstractString)
        algorithm = Symbol(algorithm)
    end

    assignment = if algorithm == :louvaincpp
        modularity_cluster(SNN.array, algorithm=1, resolution=resolution, nrandomstarts=nrandomstarts, niterations=niterations)
    elseif algorithm == :louvain
        n = Network(SNN.array)
        cl = louvain_multi(n, nrandomstarts, resolution=resolution, max_iterations=niterations)
        cl.nodecluster
    else
        error("unknown clustering algorithm: $algorithm")
    end

    NamedArray(assignment, (SNN.dicts[1],), (:cells,))
end

