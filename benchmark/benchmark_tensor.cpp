//
// Created by Roman Ellerbrock on 12/15/21.
//

#include "benchmark_tensor.h"
#include "Tensor/Tensor"
#include "Benchmark.h"

namespace benchmarks {

	Tensorcd create(size_t dim, size_t order) {
		vector<size_t> dims(order, dim);
		TensorShape shape(dims);
		return arangecd(shape);
	}

	void moveReturn() {

	}

	void swapRef() {
		size_t n = 10;
		Tensorcd mat = create(n, 2);
		Tensorcd A = create(n, 8);

		Tensorcd Tmp(A.shape_);
		for (size_t i = 0; i < A.shape_.order(); ++i) {
			matrixTensor(Tmp, mat, A, i);
			swap(A, Tmp);
		}
	}

	void swapRV() {
		size_t n = 10;
		Tensorcd mat = create(n, 2);
		Tensorcd A = create(n, 8);

		Tensorcd Tmp(A.shape_);
		cout << "critical --->\n";
		for (size_t i = 0; i < A.shape_.order(); ++i) {
			Tmp = matrixTensor(mat, A, i);
			swap(A, Tmp);
		}
		cout << "<---\n";

	}

	void runTensor() {
		benchmark(swapRef, "move semantics", 1);
		benchmark(swapRV, "swap preallocated mem", 1);
	}

}