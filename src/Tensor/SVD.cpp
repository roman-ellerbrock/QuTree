//
// Created by Roman Ellerbrock on 11/18/21.
//

#include "Tensor/SVD.h"

template <typename T>
SVD<T>::SVD(const TensorShape& shape) {
	size_t a = shape.lastBefore();
	size_t b = shape.lastDimension();
	assert(a >= b);
	assert(shape.order() == 2);
	size_t min_val = (a < b) ? a : b;

	TensorShape left({a, min_val});
	TensorShape right({min_val, b});
	TensorShape sv({min_val});

	*this = SVD<T>({Tensor<T>(left), Tensor<T>(right), Tensord(sv)});
}

template class SVD<complex<double>>;
template class SVD<double>;

