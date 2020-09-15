#include <cstdlib>
#include <cmath>

#include <ANN/ANN.h>

extern "C" void FindNeighbours(double *data, int n, int d, int k, double eps, int *nn_index, double *distances) {
	ANNpointArray pts = annAllocPts(n,d);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < d; ++j) {
			pts[i][j] = data[i + j*n];
		}
	}

	ANNkd_tree tree(pts, n, d);
	ANNidxArray nn_idx = new ANNidx[k];
	ANNdistArray dists = new ANNdist[k];

	for(int i = 0; i < n; i++) {
		tree.annkSearch(pts[i], k, nn_idx, dists, eps);

		for (int j = 0; j < k; j++) {
			int ptr = i + j*n;
			distances[ptr] = ANN_ROOT(dists[j]);
			nn_index[ptr] = nn_idx[j] + 1;
		}
	}

	delete [] dists;
	delete [] nn_idx;
	annDeallocPts(pts);
}
