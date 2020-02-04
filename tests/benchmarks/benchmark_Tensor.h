//
// Created by Roman Ellerbrock on 2/3/20.
//

#ifndef BENCHMARK_TENSOR_H
#define BENCHMARK_TENSOR_H
#include "Core/Tensor.h"
#include <random>
#include "Core/Tensor_Extension.h"

namespace benchmark_tensor {
	TensorDim make_TensorDim(size_t order, size_t dim);

	void dot_product(mt19937& gen, size_t dim = 5, size_t max_order = 5, ostream& os = cout);

	void run();


}
#endif //BENCHMARK_TENSOR_H
