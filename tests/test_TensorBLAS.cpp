//
// Created by Roman Ellerbrock on 5/25/21.
//

#include "Core/TensorBLAS.h"
#include <UnitTest++/UnitTest++.h>
#include "Core/Tensor_Extension.h"


using namespace std;

SUITE (Tensor) {
	double eps = 1e-7;

	class Factory {
	public:
		Factory()
			: gen_(129123), dim_(16) {
			S_ = Create(dim_);
		}

		Matrixcd Create(size_t dim) {
			Matrixcd x(dim, dim);
			Tensor_Extension::generate(x, gen_);
			return x;
		}

		size_t dim_;
		Matrixcd S_;
		mt19937 gen_;

	};

	TEST_FIXTURE(Factory, Transpose) {
		Matrixcd ST1(dim_, dim_);
		Matrixcd ST2(dim_, dim_);
		transpose(&ST1[0], &S_[0], dim_, dim_);
		size_t blocksize = 4;
		transpose2(&ST2[0], &S_[0], dim_, dim_, blocksize);
			CHECK_CLOSE(0., residual(ST1, ST2), eps);

		/// === non even dim ===
		size_t d = 23;
		auto X = Create(d);
		Matrixcd XT1(d, d);
		Matrixcd XT2(d, d);
		transpose(&XT1[0], &X[0], d, d);
		transpose2(&XT2[0], &X[0], d, d, blocksize);
			CHECK_CLOSE(0., residual(XT1, XT2), eps);
	}
}

