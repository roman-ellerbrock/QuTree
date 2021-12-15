//
// Created by Roman Ellerbrock on 12/15/21.
//

#include "benchmark_tensor.h"
#include "Tensor/Tensor"

void moveInPlace() {
	Tensorcd mat({2, 2});
	Tensorcd A({2, 2, 2});

	for (size_t i = 0; i < A.shape_.order(); ++i) {
		A = matrixTensor(mat, A, i);
	}

	Tensorcd Tmp({2, 2, 2});
	Tensorcd *a = &A;
	Tensorcd *tmp = &Tmp;

	for (size_t i = 0; i < A.shape_.order(); ++i) {
		matrixTensor(*tmp, mat, *a, i);
		swap(a, tmp);
	}

}