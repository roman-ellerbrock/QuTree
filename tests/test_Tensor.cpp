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
			CreateTensors();
		}

		Tensorcd A;
		Tensorcd B;
		Tensorcd C_;
		Tensorcd C2_;
		TensorShape shape_c_;

		void CreateTensorA() {
			TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
			A = Tensorcd(tdim);
			for (size_t i = 0; i < tdim.totalDimension(); ++i) {
				A(i) = i;
			}
		}

		void CreateTensorB() {
			TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
			B = Tensorcd(tdim);
			for (size_t i = 0; i < tdim.totalDimension(); ++i) {
				B(i) = i % 3;
			}
		}

		void CreateTensorC() {
			shape_c_ = TensorShape({2, 2, 2});
			// C
			C_ = Tensorcd(shape_c_);
			for (size_t bef = 0; bef < shape_c_.before(1); ++bef) {
				for (size_t act = 0; act < shape_c_[1]; ++act) {
					for (size_t aft = 0; aft < shape_c_.after(1); ++aft) {
						C_(bef, act, aft, 1) = (double) act * QM::im;
					}
				}
			}
			// C2
			C2_ = Tensorcd(shape_c_);
			for (size_t bef = 0; bef < shape_c_.before(1); ++bef) {
				for (size_t act = 0; act < shape_c_[1]; ++act) {
					for (size_t aft = 0; aft < shape_c_.after(1); ++aft) {
						C2_(bef, act, aft, 1) = (double) aft * QM::im;
					}
				}
			}
		}

		void CreateTensors() {
			CreateTensorA();
			CreateTensorB();
			CreateTensorC();
		}
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

}
