//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "Operator/LeafMatrix.h"
#include "Tree/Tree.h"
#include "Tensor/Tensor"

Tensorcd toMatrix(const LeafOperatorcd& h, const Leaf& leaf) {
	size_t mode = leaf.par().mode_;
	size_t dim = leaf.par().dim_;
	TensorShape shape{dim, dim};
	Tensorcd v(shape);
	Tensorcd hv(shape);
	for (size_t i = 0; i < dim; ++i) {
		v(i, i) = 1.;
	}
	h.apply(leaf.basis_, hv, v);
	return contraction(v, hv);
}

Tensord toMatrix(const LeafOperator<double>& h, const Leaf& leaf) {
	size_t mode = leaf.par().mode_;
	size_t dim = leaf.par().dim_;
	TensorShape shape{dim, dim};
	Tensord v(shape);
	Tensord hv(shape);
	for (size_t i = 0; i < dim; ++i) {
		v(i, i) = 1.;
	}
	h.apply(leaf.basis_, hv, v);
	return contraction(v, hv);
}

template<typename T>
void LeafMatrix<T>::apply(const BasisAPI& basis, Tensor<T>& hA,
	const Tensor<T>& A) const {
	matrixTensor(hA, h_, A, 0);
}

template<typename T>
LeafMatrix<T>::LeafMatrix(Tensor<T> h, bool adj)
: h_(move(h)) {
	if (adj) {
		h_ = adjoint(h_);
	}
}

template
class LeafMatrix<complex<double>>;

template
class LeafMatrix<double>;

