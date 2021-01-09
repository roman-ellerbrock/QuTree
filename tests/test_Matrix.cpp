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

	Matrixcd Create() {
		Matrixcd tmp(3, 3);
		for (size_t i = 0; i < tmp.Dim1(); ++i) {
			for (size_t j = 0; j < tmp.Dim2(); ++j) {
				tmp(j, i) = i + j;
			}
		}
		return tmp;
	}

	constexpr double eps = 1e-7;

	TEST_FIXTURE (MatrixFactory, Matrix_FileIO) {
		/// Test Matrix I/O
		A.Write("matrix1.tmp.dat");
		Matrixcd N("matrix1.tmp.dat");
		bool success = A == N;
			CHECK_EQUAL(success, true);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Add) {
		auto S = A + B;
		S.Write("matrix_add.dat");
		Matrixcd S_read("matrix_add.dat");
		double r = Residual(S, S);
			CHECK_CLOSE(r, 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Subst) {
		auto D = A - B;
		D.Write("matrix_subst.dat");
		Matrixcd D_read("matrix_subst.dat");
			CHECK_CLOSE(Residual(D, D_read), 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_Prod) {
		auto D = A * B;
		D.Write("matrix_prod.dat");
		Matrixcd D_read("matrix_prod.dat");
			CHECK_CLOSE(Residual(D, D_read), 0., eps);
	}


	TEST_FIXTURE (MatrixFactory, Matrix_Diagonalization) {
		auto x = A.cDiag();
		const Matrixcd& Ua = x.first;
		const Vectord& la = x.second;
		Ua.Write("matrix_cdiag_trafo.dat");
		la.Write("matrix_cdiag_ev.dat");

		Matrixcd U("matrix_cdiag_trafo.dat");
		Vectord lambda("matrix_cdiag_ev.dat");
		auto residual_U = Residual(U, x.first);
		auto residual_L = Residual(lambda, x.second);
			CHECK_CLOSE(residual_U, 0., eps);
			CHECK_CLOSE(residual_L, 0., eps);
	}

	TEST_FIXTURE (MatrixFactory, Matrix_RoF) {
		{
			// Copy asignment operator
			auto Aca = A;
			double r = Residual(A, Aca);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Copy constructor
			auto Acc(A);
			double r = Residual(A, Acc);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Move asignment operator
			auto Ama = move(Create());
			double r = Residual(A, Ama);
				CHECK_CLOSE(r, 0., eps);
		}

		{
			// Move constructor
			auto Amc(move(Create()));
				CHECK_CLOSE(Residual(A, Amc), 0., eps);
		}
	}

	TEST (Matrix_Rebuild) {
		mt19937 gen(1990);
		size_t dim = 10;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		auto x = Diagonalize(A);
		Matrixcd B = toMatrix(x);
		auto residual = Residual(A, B);
		CHECK_CLOSE(0., residual, eps);
	}

	TEST(Matrix_BuildInverse) {
		mt19937 gen(1990);
		size_t dim = 10;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		Matrixcd Adag = A.Adjoint();
		A = A * Adag;
		auto x = Diagonalize(A);

		Matrixcd B = BuildInverse(x, 1e-10);
		Matrixcd I_test = A * B;
		Matrixcd I = IdentityMatrix<complex<double>>(A.Dim2());
		auto residual = Residual(I, I_test);
		CHECK_CLOSE(0., residual, eps);
	}

	TEST(svd) {
		mt19937 gen(1293123);
		size_t dim = 20;
		Matrixcd A = RandomMatrices::GUE(dim, gen);
		auto Asvd = svd(A);
		auto B = toMatrix(Asvd);
		CHECK_CLOSE(0., Residual(A, B), 1e-8);
	}

	TEST(Toeph) {
		Matrixcd A(5, 5);
		Matrixcd B(5, 5);
		double c = 2.;

		auto C = (c * A) + B * (A + 2. * B);
		C = A + B;
		A *= c;
		B += A;
		Vectorcd v = A;
	}
}
