#include "LeafFunction.h"


template<typename T>
void LeafFunction<T>::Apply(const PrimitiveBasis& grid,
	Tensor<T>& hAcoeff, const Tensor<T>& Acoeff) const {
	h_(grid, hAcoeff, Acoeff);
}

template class LeafFunction<complex<double>>;
template class LeafFunction<double>;
