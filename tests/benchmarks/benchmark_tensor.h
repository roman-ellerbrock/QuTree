//
// Created by Roman Ellerbrock on 2/3/20.
//

#ifndef BENCHMARK_TENSOR_H
#define BENCHMARK_TENSOR_H
#include "Tensor/Tensor.h"
#include <random>
#include "Tensor/TensorLapack.h"
#include "benchmark_helper.h"

namespace benchmark {
	TensorShape make_TensorDim(size_t order, size_t dim);

	void screen_dim_contraction(mt19937& gen, ostream& os, size_t nsample, size_t mode, bool gemm_only);
	void screen_dim_matrixTensor(mt19937& gen, ostream& os, size_t nsample, size_t mode, bool gemm_only);
	void screen_order_contraction(mt19937& gen, ostream& os, size_t nsample, size_t which, bool gemm_only);
	void screen_order_matrixTensor(mt19937& gen, ostream& os, size_t nsample, size_t which, bool gemm_only);
	void benchmark_peak_gemm(mt19937& gen, ostream& os, size_t nsample);

	void run();


}
#endif //BENCHMARK_TENSOR_H
