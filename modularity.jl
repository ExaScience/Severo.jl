import SparseArrays: SparseMatrixCSC, rowvals, nzrange, nonzeros, nnz
import Random: shuffle

struct Edge
	node::Int64
	weight::Float64
end

struct Node
	self::Float64 #REMOVE
	weight::Float64
	edges::UnitRange{Int64} #REMOVE?
end

struct Network
	nodes::Vector{Node}
	edges::Vector{Edge}
	totw::Float64
end

struct Cluster
	w_in::Float64 #REMOVE modularity shouldn't be called often
	w_tot::Float64
end

mutable struct Clustering
	network::Network
	nodecluster::Vector{Int64}
	clusters::Vector{Cluster}
	nclusters::Int64
end

nodes(network::Network) = network.nodes
numnodes(network::Network) = length(network.nodes)
numedges(network::Network) = length(network.edges)
self_weight(node::Node) = node.self
numclusters(clustering::Clustering) = clustering.nclusters
numnodes(clustering::Clustering) = numnodes(clustering.network)

alledges(network::Network, nodeid::Int64) = alledges(network, network.nodes[nodeid])
alledges(network::Network, node::Node) = @inbounds view(network.edges, node.edges)

function Clustering(network::Network)
	nodecluster = collect(1:numnodes(network))
	clusters = map(nodes(network)) do node
		w_in = self_weight(node)
		w_tot = total_weight(network, node)
		Cluster(w_in, w_tot)
	end

	Clustering(network, nodecluster, clusters, numnodes(network))
end

function Clustering(network::Network, nodecluster::Vector{Int64})
	num_clusters = length(unique(nodecluster))
	clusters = map(1:num_clusters) do ci
		w_in = w_out = 0.0
		@inbounds for i in 1:numnodes(network)
			nodecluster[i] == ci || continue

			node = network.nodes[i]
			w_in += self_weight(node)
			for e in alledges(network, node)
				if nodecluster[e.node] == ci
					w_in += e.weight
				else
					w_out += e.weight
				end
			end
		end
		Cluster(w_in, w_in + w_out)
	end

	Clustering(network, nodecluster, clusters, num_clusters)
end

total_weight(network::Network, nodeid::Int64) = total_weight(network, network.nodes[nodeid])
total_weight(network::Network, node::Node) = node.weight
total_weight(network::Network) = network.totw

function cluster_weights!(kin::Vector{Float64}, neighbourcls::Vector{Int64}, clustering::Clustering, nodeid::Int64)
	# clear old
	@inbounds for neighbour in neighbourcls
		kin[neighbour] = 0.0
	end
	empty!(neighbourcls)

	network = clustering.network
	@inbounds node = network.nodes[nodeid]

	@inbounds for e in alledges(network, node)
		cj = clustering.nodecluster[e.node]
		if kin[cj] == 0.0
			push!(neighbourcls, cj)
		end

		kin[cj] += e.weight
	end

	length(neighbourcls)
end

function cluster_weights(clustering::Clustering, nodeid::Int64)
	counts = zeros(Float64, numclusters(clustering))

	network = clustering.network
	node = network.nodes[nodeid]

	@inbounds for e in alledges(network, node)
		cj = clustering.nodecluster[e.node]
		counts[cj] += e.weight
	end

	counts
end

function modularity_gain(clustering::Clustering, nodeid::Int64, to::Int64)
	totw = total_weight(clustering.network)
	kin = cluster_weights(clustering, nodeid)
	ki = total_weight(clustering.network, nodeid)

	@inbounds from = clustering.nodecluster[nodeid]
	if from == to
		0.0
	else
		@inbounds delta = (-kin[from] + (clustering.clusters[from].w_tot*ki)/totw) +
					 (kin[to] - (clustering.clusters[to].w_tot*ki)/totw) - ki^2/totw
		2delta / totw
	end
end

function modularity(clustering::Clustering)
	totw = total_weight(clustering.network)
	@inline modularity(c::Cluster) = c.w_in/totw - (c.w_tot/totw)^2
	sum(i -> modularity(clustering.clusters[i]), 1:numclusters(clustering))
end

function adjust_cluster(cluster::Cluster, kin::Float64, ki::Float64)
	w_in = cluster.w_in + 2kin
	w_out = cluster.w_tot + ki
	Cluster(w_in, w_out)
end

function move_node!(clustering::Clustering, nodeid::Int64, to::Int64)
	kin = cluster_weights(clustering, nodeid)
	ki = total_weight(clustering.network, nodeid)
	move_node!(clustering, nodeid, to, kin, ki)
end

function move_node!(clustering::Clustering, nodeid::Int64, to::Int64, kin::Vector{Float64}, ki::Float64)
	from = clustering.nodecluster[nodeid]
	clustering.clusters[from] = adjust_cluster(clustering.clusters[from], -kin[from], -ki)
	clustering.clusters[to] = adjust_cluster(clustering.clusters[to], kin[to], ki)
	clustering.nodecluster[nodeid] = to
	clustering
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
	edges_so_far = 0
	@inbounds for i in 1:nnodes
		selfweight = 0.0
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

		nodes[i] = Node(selfweight, weight, edgesstart:edges_so_far)
	end

	Network(nodes, edges, totw)
end

function reduced_network(clustering::Clustering)
	network = clustering.network
	nnodes = numclusters(clustering)
	nodes = Vector{Node}(undef, nnodes)

	edges = Vector{Edge}()
	sizehint!(edges, min(nnodes^2, numedges(network)))

	totw = 0.0
	reducedEdges = Vector{Float64}(undef, nnodes)
	for i in 1:nnodes
		fill!(reducedEdges, 0.0)

		self_weight = 0.0
		weight = 0.0

		for (nodeid, clus) in enumerate(clustering.nodecluster) #XXX
			clus == i || continue

			self_weight += network.nodes[nodeid].self
			@inbounds for e in alledges(network, nodeid)
				clust_to = clustering.nodecluster[e.node]
				if clust_to == i
					self_weight += e.weight
				else
					reducedEdges[clust_to] += e.weight
					weight += e.weight
				end
			end
		end

		edge_start = length(edges) + 1
		for (j,w) in enumerate(reducedEdges)
			w != 0.0 || continue
			push!(edges, Edge(j, w))
		end

		#XXX self_weight = w_in from cluster
		weight += self_weight
		nodes[i] = Node(self_weight, weight, edge_start:length(edges))
		totw += weight
	end

	#XXX totw = network.totw
	Network(nodes, edges, totw)
end

function Base.findmax(f, itr)
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

function best_local_move(clustering::Clustering, nodeid::Int64, neighbourcls::Vector{Int64}, kin::Vector{Float64}, ki::Float64, totw::Float64)
	@inbounds from = clustering.nodecluster[nodeid]

	(delta, idx) = findmax(neighbourcls) do neighbour_cluster
		@inbounds delta = kin[neighbour_cluster] - (clustering.clusters[neighbour_cluster].w_tot*ki)/totw
		if from == neighbour_cluster # XXX remove this
			delta += ki^2/totw
		end

		if nodeid == 0
			println("$neighbour_cluster $delta $(kin[neighbour_cluster]) $(clustering.clusters[neighbour_cluster].w_tot) $ki")
		end

		delta
	end

	@inbounds delta += -kin[from] + (clustering.clusters[from].w_tot*ki)/totw - ki^2/totw
	(2delta / totw, idx)
end

function gainz(clustering::Clustering, nodeid::Int64, neighbourcls::Vector{Int64}, kin::Vector{Float64}, ki::Float64, totw::Float64)
	@inbounds from = clustering.nodecluster[nodeid]

	delta = map(neighbourcls) do neighbour_cluster
		@inbounds kin[neighbour_cluster] - (clustering.clusters[neighbour_cluster].w_tot*ki)/totw
	end

	@inbounds delta .+= -kin[from] + (clustering.clusters[from].w_tot*ki)/totw - ki^2/totw
	2 .* delta ./ totw
end

function local_move!(clustering::Clustering; min_modularity=0.0001)
	network = clustering.network
	totw = total_weight(clustering.network)

	order = shuffle(1:numnodes(network))
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
			cluster_weights!(kin, neighbourcls, clustering, nodeid)
			ki = total_weight(clustering.network, nodeid)

			gain, bestcl = best_local_move(clustering, nodeid, neighbourcls, kin, ki, totw)
			if gain > 0.0
				move_node!(clustering, nodeid, bestcl, kin, ki)
				delta_mod += gain
				stable = false
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

	for i in 1:length(clustering.clusters)
		c, p = clustering.clusters[i], labels[i]
		while p != i && p != 0
			c, clustering.clusters[p] = clustering.clusters[p], c
			p = labels[p]
		end
	end

	for i in id:length(clustering.clusters)
		clustering.clusters[i] = Cluster(0.0, 0.0)
	end

	clustering.nclusters = id - 1
	clustering
end

function merge!(clustering::Clustering, cluster_reduced::Clustering)
	for i in 1:length(clustering.nodecluster)
		clustering.nodecluster[i] = cluster_reduced.nodecluster[clustering.nodecluster[i]]
	end
	clustering.clusters = cluster_reduced.clusters
	clustering.nclusters = cluster_reduced.nclusters
	clustering
end

function louvain!(clustering::Clustering; min_modularity=0.0001)
	numnodes(clustering) > 1 || return clustering

	gain = local_move!(clustering; min_modularity=min_modularity)
	renumber!(clustering)

	if numclusters(clustering) < numnodes(clustering)
		reduced_clustering = Clustering(reduced_network(clustering))
		gain += louvain!(reduced_clustering; min_modularity=min_modularity)
		clustering = merge!(clustering, reduced_clustering)
	end

	gain
end

function louvain(network::Network; min_modularity=0.0001, max_iterations=10)
	clustering = Clustering(network)
	mod_pre = modularity(clustering)

	total_gain = gain = louvain!(clustering; min_modularity=min_modularity)

	iter = 1
	while gain >= min_modularity && iter <= max_iterations
		gain = louvain!(clustering; min_modularity=min_modularity)
		total_gain += gain
		iter += 1
	end

	mod_post = modularity(clustering)
	println("$total_gain ≈ $(mod_post - mod_pre) $mod_post")

	clustering
end

#=
A = begin
	indices = [ 1,  2,  3,  4,  5,  6,  7,  8, 10, 11, 12, 13, 17, 19, 21, 31,  0,
	2,  3,  7, 13, 17, 19, 21, 30,  0,  1,  3,  7,  8,  9, 13, 27, 28,
	32,  0,  1,  2,  7, 12, 13,  0,  6, 10,  0,  6, 10, 16,  0,  4,  5,
	16,  0,  1,  2,  3,  0,  2, 30, 32, 33,  2, 33,  0,  4,  5,  0,  0,
	3,  0,  1,  2,  3, 33, 32, 33, 32, 33,  5,  6,  0,  1, 32, 33,  0,
	1, 33, 32, 33,  0,  1, 32, 33, 25, 27, 29, 32, 33, 25, 27, 31, 23,
	24, 31, 29, 33,  2, 23, 24, 33,  2, 31, 33, 23, 26, 32, 33,  1,  8,
	32, 33,  0, 24, 25, 28, 32, 33,  2,  8, 14, 15, 18, 20, 22, 23, 29,
	30, 31, 33,  8,  9, 13, 14, 15, 18, 19, 20, 22, 23, 26, 27, 28, 29,
	30, 31, 32]
	indptr = [  0,  16,  25,  35,  41,  44,  48,  52,  56,  61,  63,  66,  67,
	69,  74,  76,  78,  80,  82,  84,  87,  89,  91,  93,  98, 101,
	104, 106, 110, 113, 117, 121, 127, 139, 156]
	data = Float64[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1]
	dim = (34, 34)
	SparseMatrixCSC(dim[1], dim[2], indptr .+ 1, indices .+ 1, data)
end


n = Network(A)
cl = louvain(n)
println(modularity(cl))
=#
