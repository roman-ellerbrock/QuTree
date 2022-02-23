//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Tensor/Tensor_Implementation.h"

typedef double d;
typedef complex<double> cd;

template class Tensor<d>;
template class Tensor<cd>;

template Tensor<cd> random(const TensorShape& shape, mt19937& gen);
template Tensor<d> random(const TensorShape& shape, mt19937& gen);

template Tensor<cd> randomGen(const TensorShape& shape);
template Tensor<d> randomGen(const TensorShape& shape);

template Tensor<cd> arange(const TensorShape& shape);
template Tensor<d> arange(const TensorShape& shape);

template Tensor<cd> identity(const TensorShape& shape);
template Tensor<d> identity(const TensorShape& shape);

template Tensor<cd> delta(const TensorShape& shape);
template Tensor<d> delta(const TensorShape& shape);

