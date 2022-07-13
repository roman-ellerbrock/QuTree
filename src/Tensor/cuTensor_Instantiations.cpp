#include "Tensor/Tensor.hpp"
#include "Tensor/cuTensor.h"

template class Tensor<double, polymorphic::cuMemory>;
