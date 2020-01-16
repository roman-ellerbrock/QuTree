#pragma once

#include "Matrix.h"

//////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////
template<typename T>
Matrix<T>::Matrix()
	:Matrix(1, 1) {}

template<typename T>
Matrix<T>::Matrix(size_t dim1_, size_t dim2_)
	:dim1(dim1_), dim2(dim2_), size(dim1_ * dim2_),
	 coeffs(new T[dim1_ * dim2_]) {
	assert(dim1 > 0);
	assert(dim2 > 0);
	Zero();
}

// Copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix& old)
	:Matrix(old.dim1, old.dim2) {
	for (size_t i = 0; i < size; i++)
		coeffs[i] = old.coeffs[i];
}

// Move constructor
template<typename T>
Matrix<T>::Matrix(Matrix&& old) noexcept
	:dim1(old.dim1), dim2(old.dim2), size(old.dim1 * old.dim2),
	 coeffs(old.coeffs) {
	old.coeffs = nullptr;
}

// Copy Assignment Operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other) {
	//@TODO: copy dims??
	Matrix tmp(other);
	*this = move(tmp);
	return *this;
}

// Move Assignment Operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix&& other) noexcept {
	//@TODO: copy dims??// seems to be done, but maybe check again
	size = other.size;
	dim1 = other.dim1;
	dim2 = other.dim2;
	delete[] coeffs;
	coeffs = other.coeffs;
	other.coeffs = nullptr;
	return *this;
}

template<typename T>
Matrix<T>::~Matrix() {
	delete[] coeffs;
}

//////////////////////////////////////////////////////////////////////
// Bracket Operators
//////////////////////////////////////////////////////////////////////
template<typename T>
inline T& Matrix<T>::operator()(const size_t i, const size_t j) const {
	assert(i < dim1);
	assert(j < dim2);
	return coeffs[j * dim1 + i];
}

template<typename T>
inline T& Matrix<T>::operator()(const size_t i, const size_t j) {
	assert(i < dim1);
	assert(j < dim2);
	return coeffs[j * dim1 + i];
}

//////////////////////////////////////////////////////////////////////
// Fundamental Math operators
//////////////////////////////////////////////////////////////////////
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& B) {
	assert(Dim1() == B.Dim1());
	assert(Dim2() == B.Dim2());
	for (size_t j = 0; j < Dim2(); j++)
		for (size_t i = 0; i < Dim1(); i++)
			operator()(i, j) += B(i, j);
	return (*this);
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& B) {
	assert(Dim1() == B.Dim1());
	assert(Dim2() == B.Dim2());
	for (size_t j = 0; j < Dim2(); j++)
		for (size_t i = 0; i < Dim1(); i++)
			operator()(i, j) -= B(i, j);
	return (*this);
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(T coeff) noexcept {
	for (size_t i = 0; i < size; ++i) {
		coeffs[i] *= coeff;
	}
	return *this;
}


//////////////////////////////////////////////////////////////////////
// More Math operators
//////////////////////////////////////////////////////////////////////

template<typename T>
double Matrix<T>::FrobeniusNorm() const {
	double norm = 0;
	for (size_t i = 0; i < dim2; ++i) {
		for (size_t j = 0; j < dim1; ++j) {
			norm += pow(abs(operator()(j, i)), 2);
		}
	}
	return sqrt(norm);
}

template<typename T>
T Matrix<T>::Trace() const {
	assert(dim1 == dim2);
	T norm = 0;
	for (size_t i = 0; i < dim1; i++) {
		norm += operator()(i, i);
	}
	return norm;
}

template<typename T>
Matrix<T> Matrix<T>::Transpose() {
	Matrix B(dim1, dim2);
	for (size_t i = 0; i < dim1; i++)
		for (size_t j = 0; j < dim2; j++)
			B(i, j) = conjugate(operator()(j, i));
	return B;
}

template<typename T>
void Matrix<T>::rDiag(Matrix<double>& Transformation, Vector<double>& ev) {
	assert(dim1 == dim2);
	assert(ev.Dim() == dim1);
	Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>((double *) coeffs, dim1, dim2);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
	solver.compute(A);
	Eigen::MatrixXd vectors = solver.eigenvectors();
	Eigen::VectorXd eigenv = solver.eigenvalues();
	for (size_t i = 0; i < dim1; i++)
		ev(i) = eigenv(i);
	for (size_t i = 0; i < dim1; i++)
		for (size_t j = 0; j < dim2; j++)
			Transformation(j, i) = vectors(j, i);

	// Phase convention
	for (size_t i = 0; i < dim1; i++) {
		if (Transformation(0, i) < 0) {
			for (size_t j = 0; j < dim1; j++) {
				Transformation(j, i) *= -1;
			}
		}
	}
}

template<typename T>
SpectralDecomposition Matrix<T>::cDiag()const {
	Matrixcd trafo(dim1, dim2);
	Vectord ev(dim1);
	cDiag(trafo, ev);
	return pair<Matrixcd, Vectord> (trafo, ev);
}

template<typename T>
void Matrix<T>::cDiag(Matrix<complex<double>>& Transformation, Vector<double>& ev) const {
	assert(dim1 == dim2);
	assert(ev.Dim() == dim1);
	Eigen::MatrixXcd A = Eigen::Map<Eigen::MatrixXcd>((complex<double> *) coeffs, dim1, dim2);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
	solver.compute(A);
	Eigen::MatrixXcd vectors = solver.eigenvectors();
	Eigen::VectorXd eigenv = solver.eigenvalues();
	for (size_t i = 0; i < dim1; i++)
		ev(i) = eigenv(i);
	for (size_t i = 0; i < dim1; i++)
		for (size_t j = 0; j < dim2; j++)
			Transformation(j, i) = vectors(j, i);

	// Phase convention
	for (size_t i = 0; i < dim1; i++) {
		if (real(Transformation(0, i)) < 0) {
			for (size_t j = 0; j < dim1; j++) {
				Transformation(j, i) *= -1;
			}
		}
	}
}

//Inverts a complex matrix
template<typename T>
Matrix<complex<double> > Matrix<T>::cInv() const {
	assert(dim1 == dim2);

	//map everything on classes of Eigen
	Eigen::MatrixXcd A = Eigen::Map<Eigen::MatrixXcd>((complex<double> *) coeffs, dim1, dim2);

	//Invert the matrix
	Eigen::MatrixXcd B = A.inverse();

	//map the result back to MCTDH Matrices
	Matrix<complex<double> > Inverse(dim1, dim2);
	for (size_t i = 0; i < dim1; i++)
		for (size_t j = 0; j < dim2; j++)
			Inverse(j, i) = B(j, i);

	return Inverse;
}

//////////////////////////////////////////////////////////////////////
// I/O
//////////////////////////////////////////////////////////////////////
template<typename T>
void Matrix<T>::print(ostream& os) const {
	for (size_t i = 0; i < dim2; i++) {
		for (size_t j = 0; j < dim1; j++) {
			os << (*this)(j, i) << " ";
		}
		os << endl;
	}
	os << endl;
}

template<typename T>
void Matrix<T>::Write(ostream& os) const {
	for (size_t j = 0; j < dim2; j++)
		for (size_t i = 0; i < dim1; i++)
			os << operator()(i, j);
}

template<typename T>
void Matrix<T>::Read(istream& os) {
	for (size_t j = 0; j < dim2; j++)
		for (size_t i = 0; i < dim1; i++)
			os >> operator()(i, j);
}

template<typename T>
void Matrix<T>::Zero() {
	for (size_t i = 0; i < size; i++)
		coeffs[i] = 0;
}

//////////////////////////////////////////////////////////////////////
// Non-Member functions
//////////////////////////////////////////////////////////////////////
template<typename T>
Matrix<T> multAB(const Matrix<T>& A, const Matrix<T>& B) {
	assert(A.Dim2() == B.Dim1());
	Matrix<T> C(A.Dim1(), B.Dim2());
	for (size_t j = 0; j < B.Dim2(); j++) {
		for (size_t i = 0; i < A.Dim1(); i++) {
			for (size_t k = 0; k < A.Dim2(); k++) {
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return C;
}

template<typename T>
Matrix<T> multATB(const Matrix<double>& A, const Matrix<T>& B) {
	assert(A.Dim1() == B.Dim1());
	Matrix<T> C(A.Dim2(), B.Dim2());
	for (size_t j = 0; j < B.Dim2(); j++) {
		for (size_t i = 0; i < A.Dim2(); i++) {
			for (size_t k = 0; k < A.Dim1(); k++) {
				// C(i, j) += A(i, k)^T * B(k, j)
				C(i, j) += A(k, i) * B(k, j);
			}
		}
	}
	return C;
}

template<typename T>
Matrix<complex<double>> multATB(const Matrix<complex<double>>& A, const Matrix<T>& B) {
	assert(A.Dim2() == B.Dim1());
	Matrix<T> C(A.Dim1(), B.Dim2());
	for (size_t j = 0; j < B.Dim2(); j++) {
		for (size_t i = 0; i < A.Dim1(); i++) {
			for (size_t k = 0; k < A.Dim2(); k++) {
				C(i, j) += conj(A(k, i)) * B(k, j);
			}
		}
	}
	return C;
}

template<typename T>
Matrix<T> addAB(const Matrix<T>& A, const Matrix<T>& B) {
	assert(A.Dim1() == B.Dim1());
	assert(A.Dim2() == B.Dim2());
	Matrix<T> C(A.Dim1(), B.Dim2());
	for (size_t j = 0; j < A.Dim2(); j++)
		for (size_t i = 0; i < A.Dim1(); i++)
			C(i, j) = A(i, j) + B(i, j);
	return C;
}

template<typename T, typename U>
Matrix<T> multscalar(const U sca, const Matrix<T>& B) {
	Matrix<T> C(B.Dim1(), B.Dim2());
	for (size_t i = 0; i < B.Dim2(); i++)
		for (size_t j = 0; j < B.Dim1(); j++)
			C(j, i) = sca * B(j, i);
	return C;
}

template<typename T>
Matrix<T> substAB(const Matrix<T>& A, const Matrix<T>& B) {
	assert(A.Dim1() == B.Dim1());
	assert(A.Dim2() == B.Dim2());
	Matrix<T> C(A.Dim1(), B.Dim2());
	for (size_t j = 0; j < A.Dim2(); j++)
		for (size_t i = 0; i < A.Dim1(); i++)
			C(i, j) = A(i, j) - B(i, j);
	return C;
}

template<typename T>
Matrix<T> UnitarySimilarityTrafo(const Matrix<T>& A,
	const Matrix<T>& B) {
	// C=B^T*A*B
	assert(A.Dim1() == B.Dim1());
	assert(A.Dim2() == B.Dim2());
	assert(A.Dim1() == A.Dim2());
	Matrix<T> C(multAB(A, B));
	return multATB(B, C);
}

template<typename T, typename U>
Vector<T> multAB(const Matrix<U>& A, const Vector<T>& B) {
	assert(B.Dim() == A.Dim2());
	Vector<T> C(A.Dim1());
	for (size_t i = 0; i < A.Dim1(); i++) {
		for (size_t j = 0; j < A.Dim2(); j++) {
			C(i) += A(i, j) * B(j);
		}
	}
	return C;
}

template<typename T, typename U>
Vector<T> multATB(const Matrix<U>& A, const Vector<T>& B) {
	assert(B.Dim() == A.Dim2());
	Vector<T> C(A.Dim1());
	for (size_t i = 0; i < A.Dim1(); i++) {
		for (size_t j = 0; j < A.Dim2(); j++) {
			C(i) += conj(A(j, i)) * B(j);
		}
	}
	return C;
}

template<typename T>
Matrix<T> Merge(const Matrix<T>& A, const Matrix<T>& B,
	const Matrix<T>& AB) {
	// Merge Block matrices into one matrix
	// 	C =	(	A	AB	)
	// 		(	AB	B	)
	size_t n = A.Dim1();
	size_t m = B.Dim1();
	assert(n == AB.Dim1());
	assert(m == AB.Dim2());
	size_t N = m + n;

	Matrix<T> C(N, N);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C(j, i) = A(j, i);
		}
	}
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < m; ++j) {
			C(j + n, i + n) = B(j, i);
		}
	}
	// Off diagonals
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C(j, i + n) = AB(j, i);
			C(i + n, j) = conj(AB(j, i));
		}
	}
	return C;
}

template<typename T>
Vectord Matrix<T>::SolveSLE(const Vectord& b_) {
	// Re-organize b to Eigen-format
	Eigen::VectorXd b(dim1);
	for (size_t i = 0; i < dim1; i++)
		b(i) = b_(i);

	// Solve equations
	Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>((double *) coeffs, dim1, dim2);
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(A);
	Eigen::VectorXd x = dec.solve(b);

	// Save to QDlib-Vectord
	Vectord x_(dim2);
	for (size_t i = 0; i < dim2; i++)
		x_(i) = x(i);
	return x_;
}

template<typename T>
Matrix<T> Regularize(const Matrix<T>& A, double eps) {
	Matrix<T> B(A);
	size_t dim = min(A.Dim1(), A.Dim2());
	for (size_t i = 0; i < dim; i++) {
		B(i, i) += eps * exp(-B(i, i) / eps);
	}
	return B;
}

template<typename T>
Matrix<T> RealSymmetrize(const Matrix<T>& A) {
	assert(A.Dim1() == A.Dim2());
	Matrix<T> Asym(A.Dim1(), A.Dim2());
	for (size_t i = 0; i < A.Dim1(); ++i) {
		for (size_t j = 0; j < A.Dim2(); ++j) {
			Asym(j, i) = (A(j, i) + A(i, j)) / 2.;
		}
	}
	return Asym;
}

template<typename T>
Matrix<T> EuclideanDistance(const Matrix<T>& A) {
	auto G = multATB(A, A);
	Matrix<T> D(G.Dim1(), G.Dim2());
	for (size_t j = 0; j < G.Dim1(); ++j) {
		for (size_t i = 0; i < G.Dim2(); ++i) {
			D(i, j) = G(i, i) + G(j, j) - 2 * G(i, j);
		}
	}
	return D;
}