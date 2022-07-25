#include "TensorTree.h"
#include "cuMemory.h"


template class TensorTree<double, polymorphic::cuMemory>;
template class TensorTree<complex<double>, polymorphic::cuMemory>;
