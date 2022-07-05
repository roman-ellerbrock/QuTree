//
// Created by Roman Ellerbrock on 6/30/22.
//

#ifndef POLYMORPHICMEMORY_H
#define POLYMORPHICMEMORY_H
#include "stdafx.h"

namespace polymorphic {

	template<typename T, class D>
	T *allocate(const size_t n);

	template<typename T, class D>
	void copy(T *dest, const T *src, const size_t n);

	template<typename T, class D>
	void free(T *data);

	struct None;

	template<typename T, class Queue = None>
	class memory {
		using size_type = size_t;
	public:
		explicit memory(size_t size);
		~memory();

		memory(memory&& B) noexcept;
		memory(const memory& B);

		memory& operator=(const memory& B);
		memory& operator=(memory&& B) noexcept;

		void resize(size_type);

		[[nodiscard]] size_type size() const { return size_; }

		T *data() { return data_; }

		const T *data() const { return data_; }

	protected:
		T *data_{nullptr};
		size_t size_{0};
	};
}


#endif //POLYMORPHICMEMORY_H
