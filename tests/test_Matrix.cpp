//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "UnitTest++/UnitTest++.h"
#include "Core/Matrix.h"
#include <random>
#include "Util/RandomMatrices.h"

SUITE (Matrix) {
	class MatrixFactory {
	public:
		MatrixFactory() {
			CreateMatrices();
		}
		Matrixcd A;
		Matrixcd B;

		void CreateMatrixA() {
			A = Matrixcd(3, 3);
			for (size_t i = 0; i < A.dim1(); ++i) {
				for (size_t j = 0; j < A.dim2(); ++j) {
					A(j, i) = i + j;
				}
			}
		}

		void CreateMatrixB() {
			B = Matrixcd(3, 3);
			for (size_t i = 0; i < B.dim1(); ++i) {
				for (size_t j = 0; j < B.dim2(); ++j) {
					B(j, i) = i * j;
				}
			}
		}

		void CreateMatrices() {
			CreateMatrixA();
			CreateMatrixB();
		}
	};

	Matrixcd Create() {
		Matrixcd tmp(3, 3);
		for (size_t i = 0; i < tmp.dim1(); ++i) {
			for (size_t j = 0; j < tmp.dim2(); ++j) {
				tmp(j, i) = i + j;
			}
		}
		return tmp;
	}

	constexpr double eps = 1e-7;

	TEST_FIXTURE (MatrixFactory, Matrix_FileIO) {
		/// Test Matrix I/O
		A.write("matrix1.tmp.dat");
		Matrixcd N("matrix1.tmp.dat");
		bool success = A == N;
			CHECK_EQUAL(success, true);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Add) {
		auto S = A + B;
		S.write("matrix_add.dat");
		Matrixcd S_read("matrix_add.dat");
		double r = residual(S, S);
			CHECK_CLOSE(r, 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Subst) {
		auto D = A - B;
		D.write("matrix_subst.dat");
		Matrixcd D_read("matrix_subst.dat");
			CHECK_CLOSE(residual(D, D_read), 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Prod) {
		auto D = A * B;
		D.write("matrix_prod.dat");
		Matrixcd D_read("matrix_prod.dat");
			CHECK_CLOSE(residual(D, D_read), 0., eps);
	}


	TEST_FIXTURE (MatrixFactory, Matrix_Diagonalization) {
		auto x = A.cDiag();
		const Matrixcd& Ua = x.first;
		const Vectord& la = x.second;
		Ua.write("matrix_cdiag_trafo.dat");
		la.write("matrix_cdiag_ev.dat");

		Matrixcd U("matrix_cdiag_trafo.dat");
		Vectord lambda("matrix_cdiag_ev.dat");
		auto residual_U = residual(U, x.first);
		auto residual_L = residual(lambda, x.second);
			CHECK_CLOSE(residual_U, 0., eps);
			CHECK_CLOSE(residual_L, 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_RoF) {
		{
			// Copy asignment operator
			auto Aca = A;
			double r = residual(A, Aca);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Copy constructor
			auto Acc(A);
			double r = residual(A, Acc);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Move asignment operator
			auto Ama = move(Create());
			double r = residual(A, Ama);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Move constructor
			auto Amc(move(Create()));
				CHECK_CLOSE(residual(A, Amc), 0., eps);
		}
	}

	TEST (Matrix_Rebuild) {
		mt19937 gen(1990);
		size_t dim = 10;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		auto x = diagonalize(A);
		Matrixcd B = toMatrix(x);
		auto res = residual(A, B);
		CHECK_CLOSE(0., res, eps);
	}

	TEST(Matrix_BuildInverse) {
		mt19937 gen(1990);
		size_t dim = 10;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		Matrixcd Adag = A.adjoint();
		A = A * Adag;
		auto x = diagonalize(A);

		Matrixcd B = toMatrix(inverse(x, 1e-10));
		Matrixcd I_test = A * B;
		Matrixcd I = identityMatrix<complex<double>>(A.dim2());
		auto res = residual(I, I_test);
		CHECK_CLOSE(0., res, eps);
	}

	TEST(svd) {
		mt19937 gen(1293123);
		size_t dim = 20;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		auto Asvd = svd(A);
		auto B = toMatrix(Asvd);
		CHECK_CLOSE(0., residual(A, B), 1e-8);
	}

}
