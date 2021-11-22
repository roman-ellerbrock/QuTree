#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Util/QMConstants.h"
#include "Tensor/Tensor.h"
#include "Tensor/TensorBLAS1.h"

using namespace std;

SUITE (Tensor) {
	class TensorFactory {
	public:
		TensorFactory() {
			A = arangecd({2, 3, 4, 2});
		}

		Tensorcd A;
	};

	Tensorcd NewTensor() {
		TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
		Tensorcd tmp(tdim);
		for (size_t i = 0; i < tdim.totalDimension(); ++i) {
			tmp(i) = i;
		}
		return tmp;
	}

	constexpr double eps = 1e-7;

	TEST (Tensor_Constructor) {
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
		Tensorcd A(tdim);
			CHECK_CLOSE(3 * 4 * 5 * 2, A.shape_.totalDimension(), eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_RoF) {
		{
			// Copy constructor
			auto Acc(A);
				CHECK_CLOSE(0., residual(A, Acc), eps);
		}

		{
			// Copy asignment operator
			auto Aca = A;
				CHECK_CLOSE(0., residual(A, Aca), eps);
		}

		{
			// Move asignment operator
			auto Ama = move(NewTensor());
				CHECK_CLOSE(0., residual(A, Ama), eps);
		}

		{
			// Move constructor
			auto Amc(NewTensor());
				CHECK_CLOSE(0., residual(A, Amc), eps);
		}
	}

	TEST (TensorDim_FileIO) {
		/// Create a TensorDim, write to file, read in again
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
		tdim.write("tdim.dat");
		TensorShape odim("tdim.dat");
		bool same = odim == tdim;
			CHECK_EQUAL(same, true);
	}

	TEST (TensorDim_Getters) {
		/// Check Getters and initialization
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
			CHECK_EQUAL(3 * 4 * 5 * 2, tdim.totalDimension());
			CHECK_EQUAL(3 * 4 * 5, tdim.lastBefore());
			CHECK_EQUAL(2, tdim.lastDimension());
	}


	TEST_FIXTURE (TensorFactory, Tensor_FileIO) {
		/// Test Tensor I/O
		A.write("tensor1.dat");
		Tensorcd B("tensor1.dat");
			CHECK_CLOSE(0., residual(A, B), eps);
	}

	TEST (Identity) {
		Tensord Id = identityd({3, 3});
		Tensord Rd({3, 3});
		Rd(0, 0) = 1.;
		Rd(1, 1) = 1.;
		Rd(2, 2) = 1.;
			CHECK_CLOSE(0., residual(Rd, Id), eps);

		Tensorcd Icd = identitycd({3, 3});
		Tensorcd Rcd({3, 3});
		Rcd(0, 0) = 1.;
		Rcd(1, 1) = 1.;
		Rcd(2, 2) = 1.;
			CHECK_CLOSE(0., residual(Rcd, Icd), eps);
	}
}
