//
// Created by Roman Ellerbrock on 11/18/21.
//

#include "Tensor/SVD.h"

template <typename T>
SVD<T>::SVD(const TensorShape& shape, size_t k) {
	size_t a = shape[k];
	size_t b = shape.before(k) * shape.after(k);
	assert(a <= b);
//	size_t min_val = (a < b) ? a : b;

	const TensorShape& left(shape);
	TensorShape right({a, a});
	TensorShape sv({a});
//	TensorShape right({min_val, b});
//	TensorShape sv({min_val});

	*this = SVD<T>({Tensor<T>(left), Tensor<T>(right), Tensord(sv)});
}

template class SVD<complex<double>>;
template class SVD<double>;

