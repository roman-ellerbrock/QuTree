#pragma once
#include "Tensor/Tensor.h"
#include "cuMemory.h"

template <typename T>
using cuTensor = Tensor<T, polymorphic::cuMemory>;

using cuTensord = cuTensor<double>;
using cuTensorcd = cuTensor<complex<double>>;
