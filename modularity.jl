import SparseArrays: SparseMatrixCSC, rowvals, nzrange, nonzeros, nnz
import Random: shuffle

struct Edge
	node::Int64
	weight::Float64
end

struct Node
	self::Float64
	weight::Float64
	edges::UnitRange{Int64}
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

struct Clustering
	network::Network
	nodecluster::Vector{Int64}
	clusters::Vector{Cluster}
end

nodes(network::Network) = network.nodes
numnodes(network::Network) = length(network.nodes)
numedges(network::Network) = length(network.edges)
self_weight(node::Node) = node.self
numclusters(clustering::Clustering) = length(clustering.clusters)

function Clustering(network::Network)
	nodecluster = collect(1:numnodes(network))
	clusters = map(nodes(network)) do node
		w_in = self_weight(node)
		w_tot = total_weight(network, node)
		Cluster(w_in, w_tot)
	end

	Clustering(network, nodecluster, clusters)
end

function Clustering(network::Network, nodecluster::Vector{Int64})
	num_clusters = length(unique(nodecluster))
	clusters = map(1:num_clusters) do ci
		w_in = w_out = 0.0
		@inbounds for i in 1:numnodes(network)
			nodecluster[i] == ci || continue

			node = network.nodes[i]
			w_in += self_weight(node)
			for e in view(network.edges, node.edges)
				if nodecluster[e.node] == ci
					w_in += e.weight
				else
					w_out += e.weight
				end
			end
		end
		Cluster(w_in, w_in + w_out)
	end

	Clustering(network, nodecluster, clusters)
end

total_weight(network::Network, nodeid::Int64) = total_weight(network, network.nodes[nodeid])
total_weight(network::Network, node::Node) = node.weight
total_weight(network::Network) = network.totw

function cluster_weights(clustering::Clustering, nodeid::Int64)
	counts = zeros(Float64, numclusters(clustering))

	network = clustering.network
	node = network.nodes[nodeid]
	ci = clustering.nodecluster[nodeid]
#	counts[ci] += self_weight(node)

	@inbounds for e in view(network.edges, node.edges)
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
	sum(modularity, clustering.clusters)
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
			totw += nz
			weight += nz

			if rv == i
				#selfweight = nz
			else
				edges_so_far += 1
				edges[edges_so_far] = Edge(rv, nz)
			end
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

		for (node, clus) in zip(network.nodes, clustering.nodecluster)
			clus == i || continue

			@inbounds for e in view(network.edges, node.edges)
				clust_to = clustering.nodecluster[e.node]
				if clust_to == i
					self_weight += e.weight
				else
					reducedEdges[clust_to] += e.weight
				end

				weight += e.weight
			end
		end

		edge_start = length(edges) + 1
		for (j,w) in enumerate(reducedEdges)
			w != 0.0 || continue
			push!(edges, Edge(j, w))
		end

		nodes[i] = Node(self_weight, weight, edge_start:length(edges))
		totw += weight
	end

	Network(nodes, edges, totw)
end

Base.findmax(f, itr) = mapfoldl(x -> (f(x), x), _rf_findmax, itr)
_rf_findmax((fm, m), (fx, x)) = isless(fm, fx) ? (fx, x) : (fm, m)

function best_local_move(clustering::Clustering, nodeid::Int64)
	totw = total_weight(clustering.network)
	kin = cluster_weights(clustering, nodeid)
	ki = total_weight(clustering.network, nodeid)
	best_local_move(clustering, nodeid, kin, ki, totw)
end

function best_local_move(clustering::Clustering, nodeid::Int64, kin::Vector{Float64}, ki::Float64, totw::Float64)
	@inbounds from = clustering.nodecluster[nodeid]

	(delta, idx) = findmax(1:length(kin)) do neighbour_cluster
		@inbounds (neighbour_cluster != from && kin[neighbour_cluster] > 0.0) || return 0.0
		@inbounds kin[neighbour_cluster] - (clustering.clusters[neighbour_cluster].w_tot*ki)/totw
	end

	if delta != 0.0
		@inbounds delta += -kin[from] + (clustering.clusters[from].w_tot*ki)/totw - ki^2/totw
		(2delta / totw, idx)
	else
		(-1.0, 0)
	end
end

function local_move!(clustering::Clustering)
	network = clustering.network
	N = numnodes(network)
	totw = total_weight(clustering.network)

	order = shuffle(1:N)
	mod_pre = modularity(clustering)

	total_gain = 0.0
	stable = false
	while ! stable
		stable = true
		for nodeid in order
			kin = cluster_weights(clustering, nodeid) #PREALLOC
			ki = total_weight(clustering.network, nodeid)

			gain, bestcl = best_local_move(clustering, nodeid, kin, ki, totw)
			if gain > 0.0
				move_node!(clustering, nodeid, bestcl, kin, ki)
				total_gain += gain
				stable = false
			end
		end
	end

	mod_post = modularity(clustering)
	println("$total_gain ≈ $(mod_post - mod_pre)")
	total_gain
end

import SparseArrays: sparse
A = sparse(Float64[
	1  1  1  1  0  0  0  0  0
	1  1  1  0  0  0  0  0  0
	1  1  1  0  0  0  0  0  0
	1  0  0  1  1  0  0  0  0
	0  0  0  1  1  1  1  0  0
	0  0  0  0  1  1  0  0  0
	0  0  0  0  1  0  1  1  0
	0  0  0  0  0  0  1  1  1
	0  0  0  0  0  0  0  1  1
])

B = sparse([
	0.0       1.01857  0.316248  0.782568  0.0       0.0       0.0       0.0      0.0
	1.01857   0.0      1.04075   0.0       0.0       0.0       0.0       0.0      0.0
	0.316248  1.04075  0.0       0.0       0.0       0.0       0.0       0.0      0.0
	0.782568  0.0      0.0       0.0       0.724342  0.0       0.0       0.0      0.0
	0.0       0.0      0.0       0.724342  0.0       0.708958  0.171639  0.0      0.0
	0.0       0.0      0.0       0.0       0.708958  0.0       0.0       0.0      0.0
	0.0       0.0      0.0       0.0       0.171639  0.0       0.0       1.08373  0.0
	0.0       0.0      0.0       0.0       0.0       0.0       1.08373   0.0      1.61592
	0.0       0.0      0.0       0.0       0.0       0.0       0.0       1.61592  0.0
])

n = Network(B)
cl = Clustering(n)
cl_best = Clustering(n, [1, 1, 1, 2, 2, 2, 3, 3, 3])
println("$((modularity(cl), modularity(cl_best))) <-> (-0.12295272686262723, 0.529828)")
