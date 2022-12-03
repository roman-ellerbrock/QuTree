//
// Created by Roman Ellerbrock on 2/4/20.
//

#ifndef BENCHMARK_HELPER_H
#define BENCHMARK_HELPER_H
#include <chrono>
#include "Tree/Tree.h"

namespace benchmark {
	pair<double, double> statistic_helper(const vector<chrono::microseconds>& dur);
	pair<double, double> statistic_helper(const vector<double>& dur);

	TensorShape make_TensorDim(size_t order, size_t dim);

	template <typename T>
	double calc_gflops(const TensorShape& dim, size_t k);

	template <typename T>
	double calc_mem(const TensorShape& dim, size_t k);

	template <typename T>
	double calc_mem_matrixTensor(const TensorShape& dim, size_t k);

}

#endif //BENCHMARK_HELPER_H
