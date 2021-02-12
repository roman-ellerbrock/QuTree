//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "TreeOperators/LeafMatrix.h"

Matrixcd toMatrix(const LeafOperatorcd& h, const Leaf& leaf) {
	size_t mode = leaf.mode();
	size_t dim = leaf.dim();
	TensorShape shape{dim, dim};
	Tensorcd v(shape);
	Tensorcd hv(shape);
	for (size_t i = 0; i < dim; ++i) {
		v(i, i) = 1.;
	}
	h.apply(leaf.interface(), hv, v);
	Matrixcd mat = v.dotProduct(hv);
	return mat;
}


template<typename T>
void LeafMatrix<T>::apply(const LeafInterface& grid, Tensor<T>& hAcoeff,
	const Tensor<T>& Acoeff) const {
	matrixTensor(hAcoeff, h_, Acoeff, 0, true);
}

template<typename T>
LeafMatrix<T>::LeafMatrix(Matrix<T> h, bool adjoint)
: h_(move(h)) {
	if (adjoint) {
		h_ = h_.adjoint();
	}
}

template
class LeafMatrix<complex<double>>;

template
class LeafMatrix<double>;

