//
// Created by Roman Ellerbrock on 6/30/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tensor/TensorAllocator.h"

using namespace std;

SUITE (TensorAllocator) {

	TEST (CPU_allocate) {
		TensorAllocator<double> A(20);
	}
}
