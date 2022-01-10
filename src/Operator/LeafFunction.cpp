#include "Operator/LeafFunction.h"


template<typename T>
void LeafFunction<T>::apply(const BasisAPI& basis,
	Tensor<T>& hA, const Tensor<T>& A) const {
	h_(basis, hA, A);
}

template class LeafFunction<complex<double>>;
template class LeafFunction<double>;
