#include "Operator/ProductOperator.h"

template<typename T>
void ProductOperator<T>::apply(Tensor<T>& hA, Tensor<T> A, const Leaf& leaf) const {
	size_t mode = leaf.par().mode_;
	for (size_t i = 0; i < targetLeaves_.size(); ++i) {
		if (targetLeaves_[i] == mode) {
			leafOperators_[i]->apply(leaf.basis_, hA, A);
			std::swap(hA, A);
		}
	}
	std::swap(hA, A);
}

template<typename T>
TensorTree<T> ProductOperator<T>::apply(TensorTree<T> Psi) const {
	applyReference(Psi);
	return Psi;
}

template<typename T>
void ProductOperator<T>::applyFactor(TensorTree<T>& Psi, size_t i) const {
	size_t leafIdx = targetLeaves_[i];
	const Leaf* leaf = Psi.leaves_[leafIdx];
	const Node& node = leaf->parent();
	Tensor<T>& A = Psi[node];
	Tensor<T> hA(A.shape_);
	leafOperators_[i]->apply(leaf->basis_, hA, A);
	Psi[node] = hA;
}

template<typename T>
void ProductOperator<T>::applyReference(TensorTree<T>& Psi) const {
	for (size_t i = 0; i < size(); ++i) {
		applyFactor(Psi, i);
	}
}

template <typename T>
ProductOperator<T> operator*(const ProductOperator<T>& A,
	const ProductOperator<T>& B) {
	ProductOperator<T> M = B;

	for (size_t i = 0; i < A.size(); i++) {
		M.push_back(A[i], A.mode(i));
	}

	return M;
}

/*template <typename T>
void ProductOperator<T>::setV(const PotentialOperator& V_) {
    v_ = V_;
    hasV_ = true;
}
*/

typedef complex<double> cd;
template class ProductOperator<complex<double>>;
template ProductOperator<cd> operator*(const ProductOperator<cd>& A, const ProductOperator<cd>& B);

template class ProductOperator<double>;
