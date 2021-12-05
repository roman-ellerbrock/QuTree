//
// Created by Roman Ellerbrock on 2/4/20.
//

#ifndef BENCHMARK_TREE_H
#define BENCHMARK_TREE_H
#include <random>
#include "TensorNetwork/deprecated/MatrixTreeFunctions.h"

namespace benchmark {
	pair<double, double> holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);

	pair<double, double> factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);

	pair<double, double> sparse_factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);

	pair<double, double> sparse_holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os);
}

#endif //BENCHMARK_TREE_H
