#include "Tensor/Tensor.hpp"
#include "Tensor/cuTensor.h"

template class Tensor<float, polymorphic::cuMemory>;
template class Tensor<double, polymorphic::cuMemory>;
template class Tensor<complex<float>, polymorphic::cuMemory>;
template class Tensor<complex<double>, polymorphic::cuMemory>;
