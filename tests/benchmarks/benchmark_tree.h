//
// Created by Roman Ellerbrock on 2/4/20.
//

#ifndef BENCHMARK_TREE_H
#define BENCHMARK_TREE_H
#include <random>
#include "TensorNetwork/contractions.h"

namespace benchmark {
/*	pair<double, double> benchmark_matrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);
*/
	void screen_nleaves(mt19937& gen, ostream& os, size_t nsample,
	size_t min_order, size_t max_order, size_t stepk);

/*	pair<double, double> factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);

	pair<double, double> sparse_factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);

	pair<double, double> sparse_holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);
		*/
}

#endif //BENCHMARK_TREE_H
