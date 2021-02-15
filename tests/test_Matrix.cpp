//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "UnitTest++/UnitTest++.h"
//#include "Core/Matrix.h"
#include "Core/Matrix_Implementation.h"
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
		Matrix<int> C_;
		Matrixcd realH2x2_;
		Matrixcd realG2x2_;
		Matrixcd imagH2x2_;
		Matrixcd imagG2x2_;
		Vectorcd v_;

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

		void CreateMatrixC() {
			C_ = Matrix<int>(3, 3);
			for (size_t i = 0; i < C_.dim1(); ++i) {
				for (size_t j = 0; j < C_.dim2(); ++j) {
					C_(j, i) = j * 3 + i; /// ordered according to output
				}
			}
		}

		void CreateMatrices() {
			CreateMatrixA();
			CreateMatrixB();
			CreateMatrixC();

			realH2x2_ = Matrixcd(2, 2);
			realH2x2_(0, 0) = 1;
			realH2x2_(1, 0) = 2;
			realH2x2_(0, 1) = 2;
			realH2x2_(1, 1) = 1;

			complex<double> im(0., 1.);
			imagH2x2_ = im * realH2x2_;

			realG2x2_ = Matrixcd(2, 2);
			realG2x2_(0, 0) = 1;
			realG2x2_(1, 0) = 2;
			realG2x2_(0, 1) = 3;
			realG2x2_(1, 1) = 4;

			imagG2x2_ = im * realG2x2_;

			v_ = Vectorcd(2);
			v_(0) = 1;
			v_(1) = 2;
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
		Matrixcd A = RandomMatrices::gue(dim, gen);
		auto x = diagonalize(A);
		Matrixcd B = toMatrix(x);
		auto res = residual(A, B);
			CHECK_CLOSE(0., res, eps);
	}

	TEST (Matrix_BuildInverse) {
		mt19937 gen(1990);
		size_t dim = 10;
		Matrixcd A = RandomMatrices::gue(dim, gen);
		Matrixcd Adag = A.adjoint();
		A = A * Adag;
		auto x = diagonalize(A);

		Matrixcd B = toMatrix(inverse(x, 1e-10));
		Matrixcd I_test = A * B;
		Matrixcd I = identityMatrix<complex<double>>(A.dim2());
		auto res = residual(I, I_test);
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (MatrixFactory, add) {
		Matrixcd res = realH2x2_ + realH2x2_;
			CHECK_CLOSE(2., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(2., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, addEqual) {
		Matrixcd res = realH2x2_;
		res += realH2x2_;
			CHECK_CLOSE(2., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(2., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, substract) {
		Matrixcd res = realH2x2_ - realH2x2_;
			CHECK_CLOSE(0., abs(res(0, 0)), eps);
			CHECK_CLOSE(0., abs(res(1, 0)), eps);
			CHECK_CLOSE(0., abs(res(0, 1)), eps);
			CHECK_CLOSE(0., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, substractEqual) {
		Matrixcd res = realH2x2_;
		res -= realH2x2_;
			CHECK_CLOSE(0., abs(res(0, 0)), eps);
			CHECK_CLOSE(0., abs(res(1, 0)), eps);
			CHECK_CLOSE(0., abs(res(0, 1)), eps);
			CHECK_CLOSE(0., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, multiplyScalarleft) {
		Matrixcd res = 2. * realH2x2_;
			CHECK_CLOSE(2., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(2., abs(res(1, 1)), eps);
	}


	TEST_FIXTURE (MatrixFactory, multiplyScalarRight) {
		Matrixcd res = realH2x2_ * 2.;
			CHECK_CLOSE(2., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(2., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, multiplyScalarEqual) {
		Matrixcd res = realH2x2_;
		res *= 2;
			CHECK_CLOSE(2., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(2., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, divideScalarEqual) {
		Matrixcd res = realH2x2_;
		res /= 2.;
			CHECK_CLOSE(0.5, abs(res(0, 0)), eps);
			CHECK_CLOSE(1., abs(res(1, 0)), eps);
			CHECK_CLOSE(1., abs(res(0, 1)), eps);
			CHECK_CLOSE(0.5, abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, matrixVector) {
		auto res = realH2x2_ * v_;
			CHECK_CLOSE(5., abs(res(0)), eps);
			CHECK_CLOSE(4., abs(res(1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, matrixMatrix) {
		auto res = realH2x2_ * realH2x2_;
			CHECK_CLOSE(5., abs(res(0, 0)), eps);
			CHECK_CLOSE(4., abs(res(1, 0)), eps);
			CHECK_CLOSE(4., abs(res(0, 1)), eps);
			CHECK_CLOSE(5., abs(res(1, 1)), eps);
	}

	TEST_FIXTURE (MatrixFactory, equal) {
		auto Acopy = A;
			CHECK_EQUAL(true, Acopy == A);
			CHECK_EQUAL(false, A == B);
	}

	TEST_FIXTURE (MatrixFactory, notEqual) {
		auto Acopy = A;
			CHECK_EQUAL(false, Acopy != A);
			CHECK_EQUAL(true, A != B);
	}

	TEST_FIXTURE (MatrixFactory, frobeniusNorm) {
		double norm = realH2x2_.frobeniusNorm();
			CHECK_CLOSE(10., norm * norm, eps);
	}

	TEST_FIXTURE (MatrixFactory, trace) {
		complex<double> trace = realH2x2_.trace();
			CHECK_CLOSE(2., real(trace), eps);
			CHECK_CLOSE(0., imag(trace), eps);
	}

	TEST_FIXTURE (MatrixFactory, adjoint) {
		Matrixcd adj = imagH2x2_.adjoint();
		auto res = residual(imagH2x2_, -1. * adj);
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (MatrixFactory, transpose) {
		Matrixcd tra = realG2x2_.transpose();
			CHECK_CLOSE(0., abs(tra(0, 0) - 1.), eps);
			CHECK_CLOSE(0., abs(tra(1, 0) - 3.), eps);
			CHECK_CLOSE(0., abs(tra(0, 1) - 2.), eps);
			CHECK_CLOSE(0., abs(tra(1, 1) - 4.), eps);
	}

	TEST_FIXTURE (MatrixFactory, zero) {
		Matrixcd zero = realG2x2_;
		zero.zero();
			CHECK_EQUAL(0., zero(0, 0));
			CHECK_EQUAL(0., zero(1, 0));
			CHECK_EQUAL(0., zero(0, 1));
			CHECK_EQUAL(0., zero(1, 1));
	}

	TEST_FIXTURE (MatrixFactory, identityMatrix) {
		Matrixcd I = identityMatrixcd(2);
			CHECK_EQUAL(1., I(0, 0));
			CHECK_EQUAL(0., I(1, 0));
			CHECK_EQUAL(0., I(0, 1));
			CHECK_EQUAL(1., I(1, 1));
	}

	TEST_FIXTURE (MatrixFactory, conjugate) {
		Matrixcd I = imagG2x2_ + realG2x2_;
		I = I.conjugate();
			CHECK_EQUAL(complex<double>(1., -1.), I(0, 0));
			CHECK_EQUAL(complex<double>(2., -2.), I(1, 0));
			CHECK_EQUAL(complex<double>(3., -3.), I(0, 1));
			CHECK_EQUAL(complex<double>(4., -4.), I(1, 1));
	}

	TEST_FIXTURE (MatrixFactory, submatrix) {
		Matrixcd sub = subMatrix(A, 2, 2);
			CHECK_EQUAL(0., sub(0, 0));
			CHECK_EQUAL(1., sub(1, 0));
			CHECK_EQUAL(1., sub(0, 1));
			CHECK_EQUAL(2., sub(1, 1));
	}

	TEST_FIXTURE (MatrixFactory, residual) {
			CHECK_CLOSE(0., residual(realG2x2_, realG2x2_), eps);
			CHECK_CLOSE(0., residual(imagG2x2_, imagG2x2_), eps);
	}

	TEST_FIXTURE (MatrixFactory, row) {
		auto row = C_.row(1);
			CHECK_EQUAL(3, row(0));
			CHECK_EQUAL(4, row(1));
			CHECK_EQUAL(5, row(2));
	}

	TEST_FIXTURE (MatrixFactory, col) {
		auto col = C_.col(1);
			CHECK_EQUAL(1, col(0));
			CHECK_EQUAL(4, col(1));
			CHECK_EQUAL(7, col(2));
	}

	TEST (svd) {
		mt19937 gen(1293123);
		size_t dim = 20;
		Matrixcd A = RandomMatrices::gue(dim, gen);
		auto Asvd = svd(A);
		auto B = toMatrix(Asvd);
			CHECK_CLOSE(0., residual(A, B), 1e-8);
	}
}
