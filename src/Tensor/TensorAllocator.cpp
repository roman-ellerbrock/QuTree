//
// Created by Roman Ellerbrock on 6/30/22.
//

#include "TensorAllocator.h"

template <typename T>
TensorAllocator<T>::TensorAllocator(size_t size) {
}

template <typename T>
TensorAllocator<T>::~TensorAllocator() {
	free(mem_);
}

template class TensorAllocator<complex<double>>;
template class TensorAllocator<double>;
