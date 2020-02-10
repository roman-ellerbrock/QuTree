//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "LeafMatrix.h"

template<typename T>
LeafMatrix<T>::LeafMatrix(FactorMatrix<T> h)
	: h_(h) {
}

template<typename T>
void LeafMatrix<T>::Apply(const LeafInterface& grid, Tensor<T>& hAcoeff,
	const Tensor<T>& Acoeff) const {
	multAB(hAcoeff, h_, Acoeff);
}

template
class LeafMatrix<complex<double>>;

template
class LeafMatrix<double>;

