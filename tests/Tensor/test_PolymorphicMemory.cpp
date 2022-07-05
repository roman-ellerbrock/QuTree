//
// Created by Roman Ellerbrock on 6/30/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tensor/PolymorphicMemory.h"

using namespace std;

SUITE (TensorAllocator) {

	TEST (CPU_allocate) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		for (size_t i = 0; i < n; ++i) {
			A.data()[i] = 1.;
		}
			CHECK_EQUAL(n, A.size());
	}

	TEST (CPU_copyConstructor) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		polymorphic::memory<double> B(A);
		for (size_t i = 0; i < n; ++i) {
			A.data()[i] = 0.;
		}
		for (size_t i = 0; i < n; ++i) {
			B.data()[i] = 1.;
		}
			CHECK_EQUAL(n, B.size());
	}

	TEST (CPU_operatorEq0) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		polymorphic::memory<double> B = A;
		for (size_t i = 0; i < n; ++i) {
			A.data()[i] = 0.;
		}
		for (size_t i = 0; i < n; ++i) {
			B.data()[i] = 1.;
		}
			CHECK_EQUAL(n, B.size());
	}

	TEST (CPU_operatorEq1) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		polymorphic::memory<double> B(n);
		B = A;
		for (size_t i = 0; i < n; ++i) {
			A.data()[i] = 0.;
		}
		for (size_t i = 0; i < n; ++i) {
			B.data()[i] = 1.;
		}
			CHECK_EQUAL(n, B.size());
	}

	TEST (CPU_operatorEq2) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		polymorphic::memory<double>& B = A;
		B = A;
		for (size_t i = 0; i < n; ++i) {
			B.data()[i] = 1.;
		}
			CHECK_EQUAL(n, B.size());
	}

	TEST (CPU_resize0) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		size_t m = 10;
		A.resize(m);
		for (size_t i = 0; i < m; ++i) {
			A.data()[i] = 1.;
		}
			CHECK_EQUAL(m, A.size());
	}

	TEST (CPU_resize1) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		size_t m = 30;
		A.resize(m);
		for (size_t i = 0; i < m; ++i) {
			A.data()[i] = 1.;
		}
			CHECK_EQUAL(m, A.size());
	}

	TEST (CPU_resize2) {
		size_t n = 20;
		polymorphic::memory<double> A(n);
		size_t m = n;
		A.resize(m);
		for (size_t i = 0; i < m; ++i) {
			A.data()[i] = 1.;
		}
			CHECK_EQUAL(m, A.size());
	}

}

