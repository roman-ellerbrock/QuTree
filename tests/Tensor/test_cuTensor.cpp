#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Util/QMConstants.h"
#include "Tensor/cuTensor.h"

using namespace std;

SUITE (cuTensor) {

	double eps = 1e-7;

	TEST (Tensor_Constructor) {
		TensorShape tdim({3, 4, 5, 2});
		cuTensord A(tdim, false);
			CHECK_EQUAL(3 * 4 * 5 * 2, A.shape_.totalDimension());
	}

}
