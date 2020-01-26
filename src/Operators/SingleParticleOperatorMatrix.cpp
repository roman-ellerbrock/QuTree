//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "SingleParticleOperatorMatrix.h"

template<typename T>
SPOM<T>::SingleParticleOperatorMatrix(FactorMatrix<T> h)
	: h_(h) {
}

template<typename T>
void SPOM<T>::Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
	const Tensor<T>& Acoeff) const {
	multAB(hAcoeff, h_, Acoeff);
}

template class SingleParticleOperatorMatrix<complex<double>>;
template class SingleParticleOperatorMatrix<double>;

