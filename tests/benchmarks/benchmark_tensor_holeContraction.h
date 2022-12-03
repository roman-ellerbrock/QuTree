#ifndef BENCHMARK_TENSOR_ORDER_H
#define BENCHMARK_TENSOR_ORDER_H
#include "Tensor/Tensor.h"
#include <random>
#include "Tensor/TensorLapack.h"
#include "benchmark_helper.h"


namespace benchmark {
	void screen_dim_contraction(mt19937& gen, ostream& os, 
	size_t nsample, size_t mode, bool gemm_only, 
	size_t min_dim, size_t max_dim, size_t step_dim);

	void screen_order_contraction(mt19937& gen, ostream& os, 
    size_t nsample, size_t which, bool gemm_only, 
    size_t mink, size_t maxk, size_t stepk);

	auto hole_product(mt19937& gen, size_t dim, size_t order, size_t mode,
		size_t nsample, ostream& os);
}

#endif // BENCHMARK_TENSOR_ORDER_H