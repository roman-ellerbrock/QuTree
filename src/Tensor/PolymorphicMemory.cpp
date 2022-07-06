//
// Created by Roman Ellerbrock on 6/30/22.
//

#include "PolymorphicMemory.h"
#include <string.h>

namespace polymorphic {
	template<typename T, class D>
	T *allocate(const size_t n) {
		return (T *) malloc(n * sizeof(T));
	}

	template<typename T, class D>
	void copy(T *dest, const T *src, const size_t n) {
		memcpy(dest, src, n * sizeof(T));
	}

	template<typename T, class D>
	void free(T *data) {
		::free(data);
	}

	template<typename T, class D>
	memory<T, D>::memory(size_t size)
		: data_(allocate<T, D>(size)), size_(size) {
	}

	template<typename T, class D>
	memory<T, D>::~memory() {
		free<T, D>(data());
	}

	template<typename T, class D>
	memory<T, D>::memory(const memory& B)
		: memory(B.size()) {
	}

	template<typename T, class D>
	memory<T, D>::memory(memory&& B) noexcept
		: data_(std::exchange(B.data_, nullptr)),
		  size_(std::exchange(B.size_, 0)) {
	}

	template<typename T, class D>
	memory<T, D>& memory<T, D>::operator=(const memory<T, D>& B) {
		if (this == &B) {
			return *this;
		} else if (B.size_ == this->size_) {
			copy<T, D>(data(), B.data(), size());
		} else {
			memory<T, D> tmp(B);
			*this = move(tmp);
		}
		return *this;
	}

	template<typename T, class D>
	memory<T, D>& memory<T, D>::operator=(memory<T, D>&& B) noexcept {
		std::swap(data_, B.data_);
		std::swap(size_, B.size_);
		return *this;
	}

	template<typename T, class D>
	void memory<T, D>::resize(size_t newSize) {
		if (size() != newSize) {
			T* newData = allocate<T, D>(newSize);
			auto cpySize = (size() < newSize) ? size() : newSize;
			copy<T, D>(newData, data(), cpySize);
			free<T, D>(data_);

			std::swap(newData, data_);
			std::swap(newSize, size_);
		}
	}

	template
	class memory<complex<double>, CPU>;

	template
	class memory<double, CPU>;
}
