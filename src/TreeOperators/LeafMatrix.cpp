//
// Created by Roman Ellerbrock on 2020-01-24.
//

#include "LeafMatrix.h"

template<typename T>
void LeafMatrix<T>::Apply(const LeafInterface& grid, Tensor<T>& hAcoeff,
	const Tensor<T>& Acoeff) const {
	multAB(hAcoeff, h_, Acoeff, 0);
}

template<typename T>
LeafMatrix<T>::LeafMatrix(Matrix<T> h)
: h_(h) {
}

template
class LeafMatrix<complex<double>>;

template
class LeafMatrix<double>;

