# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, rowvals, nzrange, nonzeros, nnz
import Random: shuffle, default_rng, AbstractRNG

struct Edge{T <: AbstractFloat}
    node::Int64
    weight::T
end

struct Node{T <: AbstractFloat}
    weight::T
    edges::UnitRange{Int64} #REMOVE?
end

struct Network{T <: AbstractFloat}
    nodes::Vector{Node{T}}
    edges::Vector{Edge{T}}
    totw::T
    tot_self::T
end

mutable struct Clustering{T <: AbstractFloat}
    network::Network{T}
    resolution::T
    nodecluster::Vector{Int64}
    w_tot::Vector{T}
    nclusters::Int64
end

nodes(network::Network) = network.nodes
numnodes(network::Network) = length(network.nodes)
numedges(network::Network) = length(network.edges)
numclusters(clustering::Clustering) = clustering.nclusters
numnodes(clustering::Clustering) = numnodes(clustering.network)

total_weight(network::Network, nodeid::Int64) = total_weight(network, network.nodes[nodeid])
total_weight(network::Network, node::Node) = node.weight
total_weight(network::Network) = network.totw

alledges(network::Network, nodeid::Int64) = alledges(network, network.nodes[nodeid])
alledges(network::Network, node::Node) = @inbounds view(network.edges, node.edges)

function Clustering(network::Network{T}, resolution::Real) where T
    nodecluster = collect(1:numnodes(network))
    clusters = map(nodes(network)) do node
        total_weight(network, node)
    end

    resolution = convert(T, resolution)
    Clustering{T}(network, resolution, nodecluster, clusters, numnodes(network))
end

function Clustering(network::Network{T}, nodecluster::Vector{Int64}, resolution::Real) where T
    num_clusters = length(unique(nodecluster))
    clusters = map(1:num_clusters) do ci
        w_tot = zero(T)
        @inbounds for i in 1:numnodes(network)
            nodecluster[i] == ci || continue

            w_tot += total_weight(network, i)
        end
        w_tot
    end

    resolution = convert(T, resolution)
    Clustering{T}(network, resolution, copy(nodecluster), clusters, num_clusters)
end

function cluster_weights!(kin::Vector{T}, neighbourcls::Vector{Int64}, clustering::Clustering{T}, nodeid::Int64, ci::Int64) where T
    # clear old
    @inbounds for neighbour in neighbourcls
        kin[neighbour] = zero(T)
    end
    empty!(neighbourcls)

    @inbounds for e in alledges(clustering.network, nodeid)
        cj = clustering.nodecluster[e.node]
        if kin[cj] == zero(T)
            push!(neighbourcls, cj)
        end

        kin[cj] += e.weight
    end
end

function cluster_weights(clustering::Clustering{T}, nodeid::Int64) where T
    counts = zeros(T, numclusters(clustering))

    network = clustering.network
    @inbounds for e in alledges(network, nodeid)
        cj = clustering.nodecluster[e.node]
        counts[cj] += e.weight
    end

    counts
end

function sum_w_in(clustering::Clustering)
    network = clustering.network
    nodecluster = clustering.nodecluster

    sum_w_in = network.tot_self
    @inbounds for nodeid in 1:numnodes(network)
        ci = nodecluster[nodeid]

        node = network.nodes[nodeid]
        for e in alledges(network, node)
            if nodecluster[e.node] == ci
                sum_w_in += e.weight
            end
        end
    end

    sum_w_in
end

function modularity(clustering::Clustering)
    totw = total_weight(clustering.network)
    sum_w_in(clustering)/totw - clustering.resolution * sum((x/totw)^2 for x in clustering.w_tot)
end

function remove_node!(clustering::Clustering{T}, nodeid::Int64, ki::T) where T
    @inbounds from = clustering.nodecluster[nodeid]
    @inbounds clustering.w_tot[from] -= ki
    @inbounds clustering.nodecluster[nodeid] = 0
    from
end

function insert_node!(clustering::Clustering{T}, nodeid::Int64, to::Int64, ki::T) where T
    @inbounds clustering.w_tot[to] += ki
    @inbounds clustering.nodecluster[nodeid] = to
end

function checksquare(A)
    m,n = size(A)
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are $(size(A))"))
    m
end

function count_selflinks(snn::SparseMatrixCSC)
    n = checksquare(snn)
    nselflinks = 0
    @inbounds for i in 1:n
        r = nzrange(snn, i)
        j = searchsortedfirst(rowvals(snn), i, first(r), last(r), Base.Forward)
        ((j > last(r)) || (rowvals(snn)[j] != i)) && continue
        nselflinks += 1
    end
    nselflinks
end

function Network(snn::SparseMatrixCSC{T,Int64}) where T
    nnodes = checksquare(snn)
    nedges = nnz(snn) - count_selflinks(snn)

    nodes = Vector{Node{T}}(undef, nnodes)
    edges = Vector{Edge{T}}(undef, nedges)

    totw = zero(T)
    tot_self = zero(T)
    edges_so_far = 0
    @inbounds for i in 1:nnodes
        weight = zero(T)
        edgesstart = edges_so_far + 1

        r = nzrange(snn, i)
        @inbounds for j in r
            rv, nz = rowvals(snn)[j], nonzeros(snn)[j]
            rv == i && continue

            totw += nz
            weight += nz
            edges_so_far += 1
            edges[edges_so_far] = Edge(rv, nz)
        end

        nodes[i] = Node{T}(weight, edgesstart:edges_so_far)
    end

    Network{T}(nodes, edges, totw, tot_self)
end

function reduced_network(clustering::Clustering{T}) where T
    network = clustering.network
    nnodes = numclusters(clustering)
    nodes = Vector{Node{T}}(undef, nnodes)

    edges = Vector{Edge{T}}()
    sizehint!(edges, min(nnodes^2, numedges(network)))

    totw = total_weight(network)
    tot_self = network.tot_self

    #ix, counts = counting_sort(clustering.nodecluster, clustering.nclusters) XXX

    reducedEdges = Vector{T}(undef, nnodes)
    for i in 1:nnodes
        fill!(reducedEdges, zero(T))

        weight = zero(T)
        for (nodeid, clus) in enumerate(clustering.nodecluster) #XXX
            clus == i || continue

            weight += total_weight(network, nodeid)
            @inbounds for e in alledges(network, nodeid)
                clust_to = clustering.nodecluster[e.node]
                if clust_to == i
                    tot_self += e.weight
                else
                    reducedEdges[clust_to] += e.weight
                end
            end
        end

        edge_start = length(edges) + 1
        for (j,w) in enumerate(reducedEdges)
            w != zero(T) || continue
            push!(edges, Edge{T}(j, w))
        end

        nodes[i] = Node{T}(weight, edge_start:length(edges))
    end

    Network{T}(nodes, edges, totw, tot_self)
end

function best_local_move(clustering::Clustering{T}, ci::Int64, neighbourcls::Vector{Int64}, kin::Vector{T}, ki::T, totw::T) where T
    isempty(neighbourcls) && return (zero(T), ci)

    @inbounds init = kin[ci] - clustering.resolution * (clustering.w_tot[ci]*ki)/totw
    (delta, idx) = findmax(neighbourcls) do neighbour_cluster
        @inbounds kin[neighbour_cluster] - clustering.resolution * (clustering.w_tot[neighbour_cluster]*ki)/totw
    end

    (2(delta - init) / totw, idx)
end

function local_move!(rng::AbstractRNG, clustering::Clustering{T}; min_modularity::T=0.0001) where T
    network = clustering.network
    totw = total_weight(clustering.network)

    order = shuffle(rng, 1:numnodes(network))
    kin = zeros(T, numclusters(clustering))

    neighbourcls = Vector{Int64}()
    sizehint!(neighbourcls, numclusters(clustering))

    delta_mod = min_modularity
    total_gain = zero(T)
    stable = false
    while ! stable && (delta_mod >= min_modularity)
        delta_mod = zero(T)
        stable = true

        for nodeid in order
            ki = total_weight(clustering.network, nodeid)
            ci = remove_node!(clustering, nodeid, ki)

            cluster_weights!(kin, neighbourcls, clustering, nodeid, ci)

            gain, bestcl = best_local_move(clustering, ci, neighbourcls, kin, ki, totw)
            if gain > zero(T)
                insert_node!(clustering, nodeid, bestcl, ki)
                delta_mod += gain
                stable = false
            else
                insert_node!(clustering, nodeid, ci, ki)
            end
        end

        total_gain += delta_mod
    end

    total_gain
end

function renumber!(clustering::Clustering)
    labels = zeros(Int64, numclusters(clustering))
    id = 1
    for (i,c) in enumerate(clustering.nodecluster)
        if labels[c] == 0
            c = labels[c] = id
            id += 1
        else
            c = labels[c]
        end

        clustering.nodecluster[i] = c
    end

    for i in 1:length(clustering.w_tot)
        c, p = clustering.w_tot[i], labels[i]
        while p != i && p != 0
            c, clustering.w_tot[p] = clustering.w_tot[p], c
            p = labels[p]
        end
    end

    for i in id:length(clustering.w_tot)
        clustering.w_tot[i] = zero(eltype(clustering.w_tot))
    end

    clustering.nclusters = id - 1
    clustering
end

function merge!(clustering::Clustering, cluster_reduced::Clustering)
    for i in 1:length(clustering.nodecluster)
        clustering.nodecluster[i] = cluster_reduced.nodecluster[clustering.nodecluster[i]]
    end
    clustering.w_tot = cluster_reduced.w_tot
    clustering.nclusters = cluster_reduced.nclusters
    clustering
end

function louvain!(rng::AbstractRNG, clustering::Clustering{T}; min_modularity::Real=0.0001) where T
    numnodes(clustering) > 1 || return zero(T)

    min_modularity = convert(T, min_modularity)
    gain = local_move!(rng, clustering; min_modularity=min_modularity)
    renumber!(clustering)

    if numclusters(clustering) < numnodes(clustering)
        reduced_clustering = Clustering(reduced_network(clustering), clustering.resolution)
        gain += louvain!(rng, reduced_clustering; min_modularity=min_modularity)
        clustering = merge!(clustering, reduced_clustering)
    end

    gain
end

louvain!(clustering::Clustering; kw...) = louvain!(default_rng(), clustering; kw...)

function louvain(rng::AbstractRNG, network::Network;
        resolution::Real=0.8, min_modularity::Real=0.0001, max_iterations::Int64=10)

    clustering = Clustering(network, resolution)
    gain = louvain!(rng, clustering; min_modularity=min_modularity)

    iter = 1
    while gain >= min_modularity && iter <= max_iterations
        gain = louvain!(rng, clustering; min_modularity=min_modularity)
        iter += 1
    end

    clustering
end

louvain(network::Network; kw...) = louvain(default_rng(), network; kw...)

function louvain_multi(rng::AbstractRNG, network::Network, nrandomstarts::Int64; verbose::Bool=false, kw...)
    m = louvain(rng, network; kw...)
    f_m = modularity(m)
    verbose && println("1: $f_m")

    for i in 2:nrandomstarts
        x = louvain(rng, network; kw...)
        f_x = modularity(x)
        verbose && println("$i: $f_x; current best $f_m")

        if isless(f_m, f_x)
            m, f_m = x, f_x
        end
    end

    m
end

louvain_multi(network::Network, nrandomstarts::Int64; kw...) = louvain_multi(default_rng(), network, nrandomstarts; kw...)
