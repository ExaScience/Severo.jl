#include <cstdlib>
#include <cmath>

#include "annoylib.h"
#include "kissrandom.h"
#include <stdexcept>

extern "C" void FindNeighbours(double *data, int n, int d, int k, int q, int *nn_index, double *distances) {
	AnnoyIndex<int32_t, double, Euclidean, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy> index(d);

	std::vector<double> c(d);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < d; ++j)
			c[j] = data[i+j*n];

		char *errormsg;
		if (!index.add_item(i, &c[0], &errormsg))
			throw std::runtime_error(errormsg);
	}

	index.build(q);

	std::vector<int32_t> nn_idx;
	std::vector<double> dists;

	nn_idx.reserve(k);
	dists.reserve(k);
	for(int i = 0; i < n; i++) {
		nn_idx.clear();
		dists.clear();

		index.get_nns_by_item(i, k, -1, &nn_idx, &dists);

		for(int j = 0; j < k; j++) {
			int ptr = i + j*n;
			distances[ptr] = dists[j];
			nn_index[ptr] = nn_idx[j] + 1;
		}
	}
}
