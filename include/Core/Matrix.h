#pragma once
#include "stdafx.h"
#include <Eigen/Dense>
#include "Vector.h"


template<typename T>
class Matrix {
/**
 * \class Matrix
 * \ingroup Core
 * \brief This class represents a Matrix.
 *
 * This is a simple Matrix-class. The functions that
 * are implemented are specially suited for QuantumDynamics
 * simulations.
 *
 * Usage:
 * Matrixcd M(dim1_, dim2_);
 * Matrixcd A(M);
 * A = M * (A + M);
 *
 * */
public:
	//////////////////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////////////////

	Matrix();

	Matrix(size_t dim1_, size_t dim2_);

	explicit Matrix(istream& is);

	explicit Matrix(const string& filename);

	// Copy constructor
	Matrix(const Matrix& old);

	// Move constructor
	Matrix(Matrix&& old)noexcept;

	// Copy Assignment Operator
	Matrix& operator=(const Matrix& other);

	// Move Assignment Operator
	Matrix& operator=(Matrix&& other)noexcept;

	~Matrix();

	//////////////////////////////////////////////////////////////////////
	/// Backet Operators
	//////////////////////////////////////////////////////////////////////

    T& operator()(const size_t i, const size_t j)const;

	T& operator()(const size_t i, const size_t j);

	inline T& operator[](const size_t idx) const {
		return coeffs_[idx];
	}

	inline T& operator[](const size_t idx) {
		return coeffs_[idx];
	}

	//////////////////////////////////////////////////////////////////////
	/// Arithmetic
	//////////////////////////////////////////////////////////////////////

	template<typename U>
	friend Matrix operator*(const U sca, const Matrix& B) {
		return multscalar(sca, B);
	}

	friend Vector<T> operator*(const Matrix& A, const Vector<T> v) {
		return multAB(A, v);
	}

	friend Matrix operator*(const Matrix& A, const Matrix& B) {
		return multAB(A, B);
	}

	friend Matrix operator+(const Matrix& A, const Matrix& B) {
		return addAB(A, B);
	}

	friend Matrix operator-(const Matrix& A, const Matrix& B) {
		return substAB(A, B);
	}

	Matrix& operator+=(const Matrix<T>& B);

	Matrix& operator-=(const Matrix<T>& B);

	Matrix& operator*=(T coeff) noexcept;
	Matrix& operator/=(T coeff) noexcept;


	bool operator==(const Matrix<T>& A)const;
	bool operator!=(const Matrix<T>& A)const;

	//////////////////////////////////////////////////////////////////////
	/// More Math operators
	//////////////////////////////////////////////////////////////////////

	// Calculate Frobenius Norm of Matrix
	double FrobeniusNorm()const;

	// Calculate the Trace
	T Trace() const;

	// Adjoint of a matrix
	Matrix Adjoint() const;

	// Transpose a matrix
	Matrix Transpose() const;

	// Diagonals of a matrix
	Vector<T> diag() const;

	// Inverts a complex matrix
	Matrix<complex<double> > cInv() const;

	// Diagonalization
	void rDiag(Matrix<double>& Transformation, Vector<double>& ev) const;

	pair<Matrix<double>, Vectord> rDiag()const;

	void cDiag(Matrix<complex<double>>& Transformation, Vector<double>& ev) const;

	pair<Matrix<complex<double>>, Vectord> cDiag()const;

	// Solve System of linear equations for Matrix<double>
	Vectord SolveSLE(const Vectord& b_);

	void Zero();

	//////////////////////////////////////////////////////////////////////
	// I/O
	//////////////////////////////////////////////////////////////////////

	void print(ostream& os = cout) const;

	void Write(const string& filename) const;

	virtual void Write(ostream& os) const;

	void Read(const string& filename);

	virtual void Read(istream& os);


	//////////////////////////////////////////////////////////////////////
	// Getter & Setter
	//////////////////////////////////////////////////////////////////////

	size_t Dim1() const { return dim1_; }

	size_t Dim2() const { return dim2_; }

	T *Coeffs() const { return coeffs_; }

protected:
	double conjugate(const double d) const {
		return d;
	}

	complex<double> conjugate(const complex<double> c) const {
		return conj(c);
	}

	T *coeffs_;
	size_t dim1_;
	size_t dim2_;
};

typedef Matrix<complex<double>> Matrixcd;
typedef Matrix<double> Matrixd;

template <typename T>
using SpectralDecomposition = pair<Matrix<T>, Vectord>;

typedef SpectralDecomposition<complex<double>> SpectralDecompositioncd;

typedef SpectralDecomposition<double> SpectralDecompositiond;

//////////////////////////////////////////////////////////////////////
// Non-Member functions
//////////////////////////////////////////////////////////////////////

template<typename T, typename U>
Vector<T> multAB(const Matrix<U>& A, const Vector<T>& B);

template<typename T, typename U>
Vector<T> multATB(const Matrix<U>& A, const Vector<T>& B);

template<typename T>
Matrix<T> multAB(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> multATB(const Matrix<double>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> Merge(const Matrix<T>& A, const Matrix<T>& B,
		const Matrix<T>& AB);

template<typename T>
Matrix<complex<double>> multATB(const Matrix<complex<double>>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> addAB(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> substAB(const Matrix<T>& A, const Matrix<T>& B);

template<typename T, typename U>
Matrix<T> multscalar(const U sca, const Matrix<T>& B);

/**
 * \brief Diagonalize a complex double matrix
 * @param A Matrix to be diagonalized
 * @return Decomposed matrix
 */
SpectralDecompositioncd Diagonalize(const Matrix<complex<double>>& A);

/**
 * \brief Diagonalize a complex double martrix
 * @param S Decomposed matrix (call-by-reference)
 * @param A Matrix to be diagonalized
 */
void Diagonalize(SpectralDecompositioncd& S,const Matrix<complex<double>>& A);

/**
 * \brief Diagonalize a double matrix
 * @param A Matrix to be diagonalized
 * @return Decomposed matrix
 */
SpectralDecompositiond Diagonalize(const Matrix<double>& A);

/**
 * \brief Diagonalize a double martrix
 * @param S Decomposed matrix (call-by-reference)
 * @param A Matrix to be diagonalized
 */
void Diagonalize(SpectralDecompositiond& S,const Matrix<double>& A);

/**
 * \brief Diagonalization routine for remaining types (Eigenvector always double; not generally applicable)
 * @tparam T
 * @param B
 * @return
 */
template<typename T>
pair<Matrix<T>, Vectord> Diagonalize(const Matrix<T>& B);

/**
 * \brief Construct matrix from its decomposition
 * @param X Decomposed matrix
 * @return Re-constructed matrix
 */
template <typename T>
Matrix<T> BuildMatrix(const SpectralDecomposition<T>& X);

/**
 * \brief Construct inverse matrix from SpectralDecomposition
 * @param X Decomposed matrix
 * @return Inverse matrix
 */
template <typename T>
Matrix<T> BuildInverse(const SpectralDecomposition<T>& X, double eps = 1e-7);

template <typename T>
SpectralDecomposition<T> inverse(SpectralDecomposition<T> X, double eps = 1e-7);

/**
 * \brief Calculate squareroot of matrix (decomposed)
 * @param X Decomposed matrix
 * @return squareroot of matrix (decomposed)
 */
template <typename T>
SpectralDecomposition<T> sqrt(SpectralDecomposition<T> X);

/**
 * \brief Create an identity matrix
 * @tparam T Type of identity matrix
 * @param dim dimension
 * @return Identity-matrix
 */
template<typename T>
Matrix<T> IdentityMatrix(size_t dim);

Matrixcd IdentityMatrixcd(size_t dim);
Matrixd IdentityMatrixd(size_t dim);

template<typename T>
Matrix<T> UnitarySimilarityTrafo(const Matrix<T>& A,
		const Matrix<T>& B);

template<typename T>
Matrix<T> EuclideanDistance(const Matrix<T>& A);

template<typename T>
double Residual(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> RealSymmetrize(const Matrix<T>& A);

template<typename T>
Matrix<T> Regularize(const Matrix<T>& A, double eps);

template<typename T>
ostream& operator<<(ostream& os, const Matrix<T>& A);
template<typename T>
istream& operator>>(istream& is, Matrix<T>& A);

Matrixcd QR(const Matrixcd& A);

Eigen::MatrixXd toEigen(Matrixd A);
Eigen::MatrixXcd toEigen(Matrixcd A);
Matrixd toQutree(Eigen::MatrixXd A);
Matrixcd toQutree(Eigen::MatrixXcd A);


