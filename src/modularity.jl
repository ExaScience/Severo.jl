import SparseArrays: SparseMatrixCSC, rowvals, nzrange, nonzeros, nnz
import Random: shuffle, default_rng, AbstractRNG

struct Edge
	node::Int64
	weight::Float64
end

struct Node
	weight::Float64
	edges::UnitRange{Int64} #REMOVE?
end

struct Network
	nodes::Vector{Node}
	edges::Vector{Edge}
	totw::Float64
	tot_self::Float64
end

mutable struct Clustering
	network::Network
	resolution::Float64
	nodecluster::Vector{Int64}
	w_tot::Vector{Float64}
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

function Clustering(network::Network, resolution::Float64)
	nodecluster = collect(1:numnodes(network))
	clusters = map(nodes(network)) do node
		total_weight(network, node)
	end

	Clustering(network, resolution, nodecluster, clusters, numnodes(network))
end

function Clustering(network::Network, nodecluster::Vector{Int64}, resolution::Float64)
	num_clusters = length(unique(nodecluster))
	clusters = map(1:num_clusters) do ci
		w_tot = 0.0
		@inbounds for i in 1:numnodes(network)
			nodecluster[i] == ci || continue

			w_tot += total_weight(network, i)
		end
		w_tot
	end

	Clustering(network, resolution, nodecluster, clusters, num_clusters)
end

function cluster_weights!(kin::Vector{Float64}, neighbourcls::Vector{Int64}, clustering::Clustering, nodeid::Int64, ci::Int64)
	# clear old
	@inbounds for neighbour in neighbourcls
		kin[neighbour] = 0.0
	end
	empty!(neighbourcls)

	@inbounds for e in alledges(clustering.network, nodeid)
		cj = clustering.nodecluster[e.node]
		if kin[cj] == 0.0
			push!(neighbourcls, cj)
		end

		kin[cj] += e.weight
	end
end

function cluster_weights(clustering::Clustering, nodeid::Int64)
	counts = zeros(Float64, numclusters(clustering))

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

function remove_node!(clustering::Clustering, nodeid::Int64, ki::Float64)
	@inbounds from = clustering.nodecluster[nodeid]
	@inbounds clustering.w_tot[from] -= ki
	@inbounds clustering.nodecluster[nodeid] = 0
	from
end

function insert_node!(clustering::Clustering, nodeid::Int64, to::Int64, ki::Float64)
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

function Network(snn::SparseMatrixCSC{Float64,Int64})
	nnodes = checksquare(snn)
	nedges = nnz(snn) - count_selflinks(snn)

	nodes = Vector{Node}(undef, nnodes)
	edges = Vector{Edge}(undef, nedges)

	totw = 0.0
	tot_self = 0.0
	edges_so_far = 0
	@inbounds for i in 1:nnodes
		weight = 0.0
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

		nodes[i] = Node(weight, edgesstart:edges_so_far)
	end

	Network(nodes, edges, totw, tot_self)
end

function reduced_network(clustering::Clustering)
	network = clustering.network
	nnodes = numclusters(clustering)
	nodes = Vector{Node}(undef, nnodes)

	edges = Vector{Edge}()
	sizehint!(edges, min(nnodes^2, numedges(network)))

	totw = total_weight(network)
	tot_self = network.tot_self

	ix, counts = counting_sort(clustering.nodecluster, clustering.nclusters)

	reducedEdges = Vector{Float64}(undef, nnodes)
	for i in 1:nnodes
		fill!(reducedEdges, 0.0)

		self_weight = 0.0
		weight = 0.0

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
			w != 0.0 || continue
			push!(edges, Edge(j, w))
		end

		nodes[i] = Node(weight, edge_start:length(edges))
	end

	Network(nodes, edges, totw, tot_self)
end

function _findmax(f, itr)
	r = iterate(itr)
	r === nothing && error("empty collection")
	m, state = r
	f_m = f(m)
	while true
		r = iterate(itr, state)
		r === nothing && break
		x, state = r
		f_x = f(x)
		if isless(f_m, f_x) || (isequal(f_m, f_x) && x < m)
			m, f_m = x, f_x
		end
	end
	(f_m, m)
end

function best_local_move(clustering::Clustering, ci::Int64, neighbourcls::Vector{Int64}, kin::Vector{Float64}, ki::Float64, totw::Float64)
	isempty(neighbourcls) && return (0.0, ci)

	@inbounds init = kin[ci] - clustering.resolution * (clustering.w_tot[ci]*ki)/totw
	(delta, idx) = _findmax(neighbourcls) do neighbour_cluster
		@inbounds kin[neighbour_cluster] - clustering.resolution * (clustering.w_tot[neighbour_cluster]*ki)/totw
	end

	(2(delta - init) / totw, idx)
end

function local_move!(rng::AbstractRNG, clustering::Clustering; min_modularity=0.0001)
	network = clustering.network
	totw = total_weight(clustering.network)

	order = shuffle(rng, 1:numnodes(network))
	kin = zeros(Float64, numclusters(clustering))

	neighbourcls = Vector{Int64}()
	sizehint!(neighbourcls, numclusters(clustering))

	delta_mod = min_modularity
	total_gain = 0.0
	stable = false
	while ! stable && (delta_mod >= min_modularity)
		delta_mod = 0.0
		stable = true

		for nodeid in order
			ki = total_weight(clustering.network, nodeid)
			ci = remove_node!(clustering, nodeid, ki)

			cluster_weights!(kin, neighbourcls, clustering, nodeid, ci)

			gain, bestcl = best_local_move(clustering, ci, neighbourcls, kin, ki, totw)
			if gain > 0.0
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
		clustering.w_tot[i] = 0.0
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

function louvain!(rng::AbstractRNG, clustering::Clustering; min_modularity=0.0001)
	numnodes(clustering) > 1 || return 0.0

	gain = local_move!(rng, clustering; min_modularity=min_modularity)
	renumber!(clustering)

	if numclusters(clustering) < numnodes(clustering)
		reduced_clustering = Clustering(reduced_network(clustering), clustering.resolution)
		gain += louvain!(reduced_clustering; min_modularity=min_modularity)
		clustering = merge!(clustering, reduced_clustering)
	end

	gain
end

louvain!(clustering::Clustering; kw...) = louvain!(default_rng(), clustering; kw...)

function louvain(rng::AbstractRNG, network::Network;
		resolution::Float64=0.8, min_modularity::Float64=0.0001, max_iterations::Int64=10)

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
