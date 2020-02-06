#include "Operators/SingleParticleOperatorFunction.h"


template<typename T>
void SingleParticleOperatorFunction<T>::Apply(const PrimitiveBasis& grid,
	Tensor<T>& hAcoeff, const Tensor<T>& Acoeff) const {
	h_(grid, hAcoeff, Acoeff);
}

template class SingleParticleOperatorFunction<complex<double>>;
template class SingleParticleOperatorFunction<double>;
