#include "Tensor.h"
#include <iostream>
#include <UnitTest++/UnitTest++.h>

using namespace std;

SUITE (Tensor) {
	class TensorFactory {
	public:
		Tensorcd A;
		Tensorcd B;

		void CreateTensorA() {
			TensorDim tdim({2, 3, 4}, 2);
			A = Tensorcd(tdim);
			for (size_t i = 0; i < tdim.getdimtot(); ++i) {
				A(i) = i;
			}
		}

		void CreateTensorB() {
			TensorDim tdim({2, 3, 4}, 2);
			B = Tensorcd(tdim);
			for (size_t i = 0; i < tdim.getdimtot(); ++i) {
				B(i) = i % 3;
			}
		}

		void CreateTensors() {
			CreateTensorA();
			CreateTensorB();
		}
	};

	constexpr double eps = 1e-7;

	TEST (TensorDim_FileIO) {
		/// Create a TensorDim, write to file, read in again
		TensorDim tdim({3, 4, 5}, 2);
		tdim.Write("tdim.dat");
		TensorDim odim("tdim.dat");
		bool same = odim == tdim;
			CHECK_EQUAL(same, true);
	}

	TEST (TensorDim_Getters) {
		/// Check Getters and Initialization
		bool success = true;
		TensorDim tdim({3, 4, 5}, 2);
		if (tdim.getdimtot() != 3 * 4 * 5 * 2) { success = false; }
		if (tdim.getdimpart() != 3 * 4 * 5) { success = false; }
		if (tdim.getntensor() != 2) { success = false; }
			CHECK_EQUAL(success, true);
	}

	TEST (Tensor_Constructor) {
		TensorDim tdim({3, 4, 5}, 3);
		Tensorcd A(tdim);
		Tensorcd B(tdim);
		auto same = A - B;
		auto s = same.DotProduct(same);
		auto delta = s.FrobeniusNorm();
			CHECK_CLOSE(delta, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_FileIO) {
		/// Test Tensor I/O
		CreateTensors();
		A.Write("tensor1.dat");
		Tensorcd B("tensor1.dat");
		Tensorcd C = A - B;
		Matrixcd s = C.DotProduct(C);
		double residual = abs(s.Trace());
			CHECK_CLOSE(residual, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Product) {
		CreateTensors();
		auto x = HoleProduct(A, B, 0);
		Matrixcd s("Tensor_Product.dat");
		auto d = x - s;
		double residual = d.FrobeniusNorm();
			CHECK_CLOSE(residual, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Matrix_Product) {
		CreateTensors();
		auto x = HoleProduct(A, B, 1);
		x.Write("Tensor_Product_0.dat");
		Matrixcd s("Tensor_Product_0.dat");
		auto d = x - s;
		double residual = d.FrobeniusNorm();
			CHECK_CLOSE(residual, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_RoF) {
		CreateTensors();
		{
			// Copy asignment operator
			auto Aca = A;
			double r = Residual(A, Aca);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Copy constructor
			auto Acc(A);
				CHECK_CLOSE(Residual(A, Acc), 0., eps);
		}

		/*
		{
			// Move asignment operator
			auto Ama = move(MoveTensor());
				CHECK_CLOSE(Residual(A, Ama), 0., eps);
		}

		{
			// Move constructor
			auto Amc(MoveTensor());
				CHECK_CLOSE(Residual(A, Amc), 0., eps);
		}
		*/
	}
}
