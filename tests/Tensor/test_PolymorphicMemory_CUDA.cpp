//
// Created by Roman Ellerbrock on 6/30/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tensor/PolymorphicMemory_CUDA.cpp"
#include <lapack.hh>

using namespace std;

SUITE (polymorphicMemoryCUDA) {

	TEST (CUDA_allocate) {
		size_t n = 20;
		polymorphic::memory<double, polymorphic::CUDA> A(n);
			CHECK_EQUAL(n, A.size());
			CHECK_EQUAL(false, (A.data() == nullptr));
	}

}

