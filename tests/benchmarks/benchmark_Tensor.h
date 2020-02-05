//
// Created by Roman Ellerbrock on 2/3/20.
//

#ifndef BENCHMARK_TENSOR_H
#define BENCHMARK_TENSOR_H
#include "Core/Tensor.h"
#include <random>
#include "Core/Tensor_Extension.h"
#include "../tests/benchmarks/benchmark_helper.h"

namespace benchmark {
	TensorDim make_TensorDim(size_t order, size_t dim);

	void run();


}
#endif //BENCHMARK_TENSOR_H
