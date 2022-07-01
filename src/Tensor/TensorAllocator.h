//
// Created by Roman Ellerbrock on 6/30/22.
//

#ifndef TENSORALLOCATOR_H
#define TENSORALLOCATOR_H
#include "stdafx.h"

template<typename T>
void copy(T *dest, const T *src);

template<typename T>
class TensorAllocator {
public:
	explicit TensorAllocator(size_t size);
	~TensorAllocator();

	TensorAllocator(TensorAllocator&& B) noexcept;
	TensorAllocator(const TensorAllocator& B);

	TensorAllocator& operator=(const TensorAllocator& B);
	TensorAllocator& operator=(TensorAllocator&& B) noexcept;

	void resize(size_t);

	T *mem() { return mem_; }
	const T *mem() const { return mem_; }

protected:
	size_t size_{0};
	T *mem_{nullptr};
	size_t used_size_{0};
};


#endif //TENSORALLOCATOR_H
