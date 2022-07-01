//
// Created by Roman Ellerbrock on 6/30/22.
//

#ifndef TENSORALLOCATOR_H
#define TENSORALLOCATOR_H
#include "stdafx.h"

template <typename T>
class TensorAllocator {
public:
	TensorAllocator(size_t size);
	~TensorAllocator();

	void resize(size_t);

protected:
	size_t size_ {0};
	T* mem_ {nullptr};
	size_t used_size_{0};
};


#endif //TENSORALLOCATOR_H
