//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "TreeOperators/LeafMatrix.h"

Matrixcd toMatrix(const LeafOperatorcd& h, const Leaf& leaf) {
	size_t mode = leaf.Mode();
	size_t dim = leaf.Dim();
	TensorShape shape{dim, dim};
	Tensorcd v(shape);
	Tensorcd hv(shape);
	for (size_t i = 0; i < dim; ++i) {
		v(i, i) = 1.;
	}
	h.Apply(leaf.PrimitiveGrid(), hv, v);
	Matrixcd mat = v.DotProduct(hv);
	return mat;
}


template<typename T>
void LeafMatrix<T>::Apply(const LeafInterface& grid, Tensor<T>& hAcoeff,
	const Tensor<T>& Acoeff) const {
	MatrixTensor(hAcoeff, h_, Acoeff, 0, true);
}

template<typename T>
LeafMatrix<T>::LeafMatrix(Matrix<T> h, bool adjoint)
: h_(move(h)) {
	if (adjoint) {
		h_ = h.Adjoint();
	}
}

template
class LeafMatrix<complex<double>>;

template
class LeafMatrix<double>;

