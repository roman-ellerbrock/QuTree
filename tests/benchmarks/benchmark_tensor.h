//
// Created by Roman Ellerbrock on 2/3/20.
//

#ifndef BENCHMARK_TENSOR_H
#define BENCHMARK_TENSOR_H
#include "Core/Tensor.h"
#include <random>
#include "Core/Tensor_Extension.h"
#include "benchmark_helper.h"

namespace benchmark {
	TensorShape make_TensorDim(size_t order, size_t dim);

	void run();


}
#endif //BENCHMARK_TENSOR_H
