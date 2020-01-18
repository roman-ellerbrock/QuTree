#include "Tensor.h"
#include <iostream>
#include <UnitTest++/UnitTest++.h>

using namespace std;

SUITE(Tensor) {
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

	TEST(TensorDim_FileIO) {
		/// Create a TensorDim, write to file, read in again
		TensorDim tdim({3, 4, 5}, 2);
		tdim.Write("tdim.dat");
		TensorDim odim("tdim.dat");
		bool same = odim == tdim;
			CHECK_EQUAL(same, true);
	}

	TEST(TensorDim_Getters) {
		/// Check Getters and Initialization
		bool success = true;
		TensorDim tdim({3, 4, 5}, 2);
		if (tdim.getdimtot() != 3 * 4 * 5 * 2) { success = false; }
		if (tdim.getdimpart() != 3 * 4 * 5) { success = false; }
		if (tdim.getntensor() != 2) { success = false; }
			CHECK_EQUAL(success, true);
	}

	TEST(Tensor_Constructor) {
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
		CreateTensorA();
		A.Write("tensor1.dat");
		Tensorcd B("tensor1.dat");
		auto C = A - B;
		auto s = C.DotProduct(C);
		auto residual = abs(s.Trace());
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
}

SUITE (Matrix) {
	class MatrixFactory {
	public:
		Matrixcd A;
		Matrixcd B;

		void CreateMatrixA() {
			A = Matrixcd(3, 3);
			for (size_t i = 0; i < A.Dim1(); ++i) {
				for (size_t j = 0; j < A.Dim2(); ++j) {
					A(j, i) = i + j;
				}
			}
		}

		void CreateMatrixB() {
			B = Matrixcd(3, 3);
			for (size_t i = 0; i < B.Dim1(); ++i) {
				for (size_t j = 0; j < B.Dim2(); ++j) {
					B(j, i) = i * j;
				}
			}
		}

		void CreateMatrices() {
			CreateMatrixA();
			CreateMatrixB();
		}
	};

	constexpr double eps = 1e-7;

	TEST_FIXTURE(MatrixFactory, Matrix_FileIO) {
		/// Test Matrix I/O
		CreateMatrices();
		A.Write("matrix1.dat");
		Matrixcd N("matrix1.dat");
		bool success = A == N;
		CHECK_EQUAL(success, true);
	}

	TEST_FIXTURE(MatrixFactory, Matrix_) {
		CreateMatrices();
		auto C = A * B;
		Matrixcd D("matrix2.dat");
		auto d = C - D;
		double residual = d.FrobeniusNorm();
		CHECK_CLOSE(residual, 0., eps);
	}
}

