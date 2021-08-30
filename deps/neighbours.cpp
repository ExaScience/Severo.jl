#include <cstdlib>
#include <cmath>

#include "annoylib.h"
#include "kissrandom.h"
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#include <shared_mutex>
#include <mutex>

class AnnoyIndexOMPBuildPolicy {
private:
  std::shared_timed_mutex nodes_mutex;
  std::mutex n_nodes_mutex;
  std::mutex roots_mutex;

public:
  template<typename S, typename T, typename D, typename Random>
  static void build(AnnoyIndex<S, T, D, Random, AnnoyIndexOMPBuildPolicy>* annoy, int q, int n_threads) {
    AnnoyIndexOMPBuildPolicy threaded_build_policy;
    if (n_threads == -1) {
      n_threads = omp_get_max_threads();
    }

	#pragma omp parallel
	{
		const int thread_idx = omp_get_thread_num();
		const int trees_per_thread = q == -1 ? -1 : (int)floor((q + thread_idx) / n_threads);
		annoy->thread_build(trees_per_thread, thread_idx, threaded_build_policy);
	}
  }

  void lock_n_nodes() {
    n_nodes_mutex.lock();
  }
  void unlock_n_nodes() {
    n_nodes_mutex.unlock();
  }

  void lock_nodes() {
    nodes_mutex.lock();
  }
  void unlock_nodes() {
    nodes_mutex.unlock();
  }

  void lock_shared_nodes() {
    nodes_mutex.lock_shared();
  }
  void unlock_shared_nodes() {
    nodes_mutex.unlock_shared();
  }

  void lock_roots() {
    roots_mutex.lock();
  }
  void unlock_roots() {
    roots_mutex.unlock();
  }
};

using BuildPolicy = AnnoyIndexOMPBuildPolicy;
#else
using BuildPolicy = AnnoyIndexSingleThreadedBuildPolicy;
#endif

struct CosineDist : Angular {
  template<typename T>
  static inline T normalized_distance(T distance) {
    return std::max(distance/2.0, T(0));
  }

  static const char* name() {
    return "cosinedist";
  }
};

template <typename D>
void FindNeighbours(double *data, int n, int d, int stride, int k, int q, bool include_self, int *nn_index, double *distances) {
	using AnnoyIndexType = AnnoyIndex<int32_t, double, D, Kiss64Random, BuildPolicy>;
	AnnoyIndexType index(d);

	std::vector<double> c(d);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < d; ++j)
			c[j] = data[i+j*stride];

		char *errormsg;
		if (!index.add_item(i, &c[0], &errormsg))
			throw std::runtime_error(errormsg);
	}

	index.build(q);

	if(!include_self)
		k = k + 1;

#ifdef _OPENMP
	#pragma omp parallel
	{
		std::vector<int32_t> nn_idx;
		std::vector<double> dists;

		nn_idx.reserve(k);
		dists.reserve(k);

		#pragma omp for
		for(int i = 0; i < n; i++) {
			nn_idx.clear();
			dists.clear();

			index.get_nns_by_item(i, k, -1, &nn_idx, &dists);

			int idx = 0;
			for(int j = 0; j < k; j++) {
				if(!include_self && nn_idx[j] == i) continue;

				int ptr = i + idx*n;
				distances[ptr] = dists[j];
				nn_index[ptr] = nn_idx[j] + 1;
				idx++;
			}
		}
	}
#else
	std::vector<int32_t> nn_idx;
	std::vector<double> dists;

	nn_idx.reserve(k);
	dists.reserve(k);

	for(int i = 0; i < n; i++) {
		nn_idx.clear();
		dists.clear();

		index.get_nns_by_item(i, k, -1, &nn_idx, &dists);

		int idx = 0;
		for(int j = 0; j < k; j++) {
			if(!include_self && nn_idx[j] == i) continue;

			int ptr = i + idx*n;
			distances[ptr] = dists[j];
			nn_index[ptr] = nn_idx[j] + 1;
			idx++;
		}
	}
#endif
}

extern "C" void FindNeighboursEuclidean(double *data, int n, int d, int stride, int k, int q, int include_self, int *nn_index, double *distances) {
	FindNeighbours<Euclidean>(data, n, d, stride, k, q, include_self != 0, nn_index, distances);
}

extern "C" void FindNeighboursCosine(double *data, int n, int d, int stride, int k, int q, int include_self, int *nn_index, double *distances) {
	FindNeighbours<CosineDist>(data, n, d, stride, k, q, include_self != 0, nn_index, distances);
}
