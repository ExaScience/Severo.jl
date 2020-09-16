#include <cstdlib>
#include <cstdio>
#include <chrono>

#include "ModularityOptimizer.h"

[[noreturn]] static void stop(const char *msg) {
	throw std::runtime_error(msg);
}

static void print_progress(int *ticks_displayed, int i, int max, FILE *out=stdout) {
	const int max_ticks = 50;
	if( i == 0 ) {
		fprintf(out, "0%%   10   20   30   40   50   60   70   80   90   100%%\n");
		fprintf(out, "[----|----|----|----|----|----|----|----|----|----|\n");
	}

	const int ticks = (double(i) / max) * max_ticks;
	const int delta = ticks - *ticks_displayed;
	if( delta > 0 ) {
		for (int i = 0; i < delta; ++i) {
			fprintf(out, "*");
			fflush(out);
		}
		*ticks_displayed = ticks;
	}

	if(*ticks_displayed >= max_ticks) {
		fprintf(out, "|\n");
	}
}

extern "C"
void ModularityClustering(int64_t *indptr, int64_t *indices, double *values, int m, int n, int64_t nnz,
		int modularityFunction, double resolution, int algorithm, int nRandomStarts, int nIterations, int randomSeed, bool verbose,
		int *assignment /* output, size max(m,n)*/) {
	using namespace ModularityOptimizer;

	if (modularityFunction != 1 && modularityFunction != 2)
		stop("Modularity parameter must be equal to 1 or 2.");
	if (algorithm != 1 && algorithm != 2 && algorithm != 3 && algorithm != 4)
		stop("Algorithm for modularity optimization must be 1, 2, 3, or 4");
	if (nRandomStarts < 1)
		stop("Have to have at least one start");
	if (nIterations < 1)
		stop("Need at least one interation");
	if (modularityFunction == 2 && resolution > 1.0)
		stop("error: resolution<1 for alternative modularity");

	if(verbose) printf("Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck\n");

	std::shared_ptr<Network> network;
	{
		// Load lower triangle
		int64_t network_size = (nnz / 2) + 3;
		IVector node1;
		IVector node2;
		DVector edgeweights;
		node1.reserve(network_size);
		node2.reserve(network_size);
		edgeweights.reserve(network_size);
		for (int col = 0; col < n; ++col) {
			for(auto idx = indptr[col]; idx < indptr[col+1]; ++idx) {
				// 1-based indexing in the sparse matrix
				int row = indices[idx - 1] - 1;
				double value = values[idx - 1];
				if(col >= row) continue;

				node1.emplace_back(col);
				node2.emplace_back(row);
				edgeweights.emplace_back(value);
			}
		}
		if (node1.size() == 0)
			stop("Matrix contained no network data.  Check format.");

		int nNodes = std::max(m, n);
		network = matrixToNetwork(node1, node2, edgeweights, modularityFunction, nNodes);
	}

	if(verbose) {
		printf("Number of nodes: %d\n", network->getNNodes());
		printf("Number of edges: %d\n", network->getNEdges());
		printf("Running %s...\n",
			((algorithm == 1) ? "Louvain algorithm" :
			((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : "smart local moving algorithm")));
	}

	using namespace std::chrono;
	double resolution2 = ((modularityFunction == 1) ? (resolution / (2 * network->getTotalEdgeWeight() + network->getTotalEdgeWeightSelfLinks())) : resolution);

	auto beginTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
	std::shared_ptr<Clustering> clustering;
	double maxModularity = -std::numeric_limits<double>::infinity();
	JavaRandom random(randomSeed);

	int progress = 0;
	for(int i = 0; i < nRandomStarts; i++) {
		VOSClusteringTechnique vosClusteringTechnique(network, resolution2);
		if(verbose) print_progress(&progress, i, nRandomStarts);

		int j = 0;
		bool update = true;
		double modularity = 0.0;
		do {
			if (algorithm == 1)
				update = vosClusteringTechnique.runLouvainAlgorithm(random);
			else if (algorithm == 2)
				update = vosClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
			else if (algorithm == 3)
				vosClusteringTechnique.runSmartLocalMovingAlgorithm(random);
			j++;

			modularity = vosClusteringTechnique.calcQualityFunction();
		} while ((j < nIterations) && update);

		if(modularity > maxModularity) {
			clustering = vosClusteringTechnique.getClustering();
			maxModularity = modularity;
		}
	}
	if(verbose) print_progress(&progress, nRandomStarts, nRandomStarts);

	auto endTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
	if (clustering == nullptr)
		stop("Clustering step failed.");

	if(verbose) {
		if (nRandomStarts == 1) {
			if (nIterations > 1)
				printf("\n");
			printf("Modularity: %.4f\n", maxModularity);
		}
		else
			printf("Maximum modularity in %d random starts: %.4f\n", nRandomStarts, maxModularity);
		printf("Number of communities: %d\n", clustering->getNClusters());
		printf("Elapsed time: %d seconds\n", static_cast<int>((endTime - beginTime).count() / 1000.0));
	}

	// Return results
	clustering->orderClustersByNNodes();
	std::copy(clustering->cluster.cbegin(), clustering->cluster.cend(), assignment);
}
