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
 * \code{.cpp}
 * Matrixcd M(dim1, dim2);
 * Matrixcd A(M);
 * A = M * (A + M);
 * \endcode
 *
 * */
public:
	//////////////////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////////////////

	/** \brief Zero constructor equivalent to Matrix(1,1);
	 */
	Matrix();

	/** \brief Zero constructor
	 *
	 * @param dim1 size of first dimension
	 * @param dim2 size of second dimension
	 */
	Matrix(size_t dim1, size_t dim2);

	/** \brief Zero constructor
	 *
	 * @param dim1 size of first dimension
	 * @param dim2 size of second dimension
	 * @param ptr pointer to pre-allocated memory
	 * @param ownership if true, matrix will take care of deallocation
	 * @param InitZero if true, set entries to zero
	 */
	Matrix(size_t dim1, size_t dim2, T *ptr, bool ownership, bool InitZero);

	/** \brief Constructor from input stream
	 *
	 * @param is Input stream passed to read()
	 */
	explicit Matrix(istream& is);

	/** \brief Constructor from file
	 *
	 * @param filename File to load and pass to read()
	 */
	explicit Matrix(const string& filename);

	/** \brief Copy constructor
	 *
	 * @param old Matrix to be copied
	 */
	Matrix(const Matrix& old);

	/** \brief Move constructor
	 *
	 * @param old Matrix to be moved
	 */
	Matrix(Matrix&& old)noexcept;

	/** \brief Copy assignment constructor
	 *
	 * @param other Matrix to copy
	 * @return Copied matrix
	 */
	Matrix& operator=(const Matrix& other);

	/** \brief Move assignment constructor
	 *
	 * @param other Matrix to move
	 * @return Moved Matrix
	 */
	Matrix& operator=(Matrix&& other)noexcept;

	/** \brief Destructor
	 */
	~Matrix();

    /////////////////////////////////////////////////////////////////////////
	/** \name Bracket Operators
	 * Group of functions for accessing Matrix elements
	 */
	///@{

	/** \brief Two-indexed accessor
	 *
	 * @param i Row
	 * @param j Column
	 * @return Value
	 */
    T& operator()(const size_t i, const size_t j)const;

    /** \brief Two-indexed accessor
     *
     * @param i Row
     * @param j Column
     * @return Value
     */
	T& operator()(const size_t i, const size_t j);

	/** \brief Compound-index accessor
	 *
	 * @param idx Compound index
	 * @return Value
	 */
	inline T& operator[](const size_t idx) const {
		return coeffs_[idx];
	}

    /** \brief Compound-index accessor
     *
     * @param idx Compound index
     * @return Value
     */
	inline T& operator[](const size_t idx) {
		return coeffs_[idx];
	}

	///@} End Bracket Operators

    /////////////////////////////////////////////////////////////////////////
    /** \name Arithmetic Operators
     * Group of functions and operator overloadings
     * for performing arithmetic operations on Matrices
     */
    ///@{

	/** \brief Scalar elementwise multiplication
	 *
	 * @tparam U Type for scalar
	 * @param sca Scalar value
	 * @param B Target Matrix
	 * @return Scaled Matrix
	 */
	template<typename U>
	friend Matrix operator*(const U sca, const Matrix& B) {
		return multscalar(sca, B);
	}
	template<typename U>
	friend Matrix operator*(const Matrix& B, const U sca) {
		return multscalar(sca, B);
	}


	/** \brief Matrix-vector product
	 *
	 * @param A Matrix
	 * @param v Vector
	 * @return Result of A*v
	 */
	friend Vector<T> operator*(const Matrix& A, const Vector<T> v) {
		return multAB(A, v);
	}

	/** \brief Matrix-matrix product
	 *
	 * @param A Left Matrix
	 * @param B Right Matrix
	 * @return Result of A*B
	 */
	friend Matrix operator*(const Matrix& A, const Matrix& B) {
		return multAB(A, B);
	}

	/** \brief Elementwise matrix addition
	 *
	 * @param A Matrix
	 * @param B Matrix
	 * @return Result of A + B
	 */
	friend Matrix operator+(const Matrix& A, const Matrix& B) {
		return addAB(A, B);
	}

	/** \brief Elementwise matrix subtraction
	 *
	 * @param A Matrix
	 * @param B Matrix
	 * @return Result of A - B
	 */
	friend Matrix operator-(const Matrix& A, const Matrix& B) {
		return substAB(A, B);
	}

	/** \brief In-place elementwise matrix addition
	 *
	 * @param B Matrix of increments
	 * @return Incremented Matrix
	 */
	Matrix& operator+=(const Matrix<T>& B);

	/** \brief In-place elementwise matrix subtraction
	 *
	 * @param B Matrix of decrements
	 * @return Decremented Matrix
	 */
	Matrix& operator-=(const Matrix<T>& B);

	/** \brief In-place elementwise scalar matrix multiplication
	 *
	 * @param coeff Scalar
	 * @return Scaled Matrix
	 */
	Matrix& operator*=(T coeff) noexcept;

    /** \brief In-place elementwise scalar matrix division
    *
    * @param coeff Scalar
    * @return Scaled Matrix
    */
	Matrix& operator/=(T coeff) noexcept;


	/** \brief Matrix equality
	 *
	 * @param A Test Matrix
	 * @return True if precisely equal, False otherwise
	 */
	bool operator==(const Matrix<T>& A)const;

	/** \brief Matrix inequality
	 *
	 * @param A Test Matrix
	 * @return True if not precisely equal, False otherwise
	 */
	bool operator!=(const Matrix<T>& A)const;

    ///@} End Arithmetic Operators

    /////////////////////////////////////////////////////////////////////////
    /** \name Math Operators
     * Group of functions for common Matrix mathematical operations
     */
    ///@{

	/**
	 * @return Frobenius norm
	 */
	double frobeniusNorm()const;

	/**
	 * @return trace of matrix
	 */
	T trace() const;

    /**
     * @return Adjoint of matrix
     */
	Matrix adjoint() const;

	/**
	 * @return conjugate of matrix
	 */
	Matrix conjugate() const;

	/**
	 * @return Transposed matrix
	 */
	Matrix transpose() const;

	/** \brief invert a complex Matrix
	 *
	 * @return Inverted Matrix
	 */
	Matrix<T> cInv() const;

	/** \brief Diagonalization of a real Matrix
	 *
	 * Uses Eigen to diagonalize the Matrix.
	 * Phases eigenvectors s.t. first element is always positive.
	 *
	 * @param Transformation Matrix of eigenvectors
	 * @param ev Vector of eigenvalues
	 */
	void rDiag(Matrix<double>& Transformation, Vector<double>& ev) const;

	/** \brief Diagonalization of a real Matrix
	 *
	 * Exactly like `Matrix::rDiag(Matrixd&, Vectord&)` but returns copies instead of pass by reference
	 *
	 * @return Eigenvectors/value as a pair of <Matrixd, Vectord>
	 */
	pair<Matrix<double>, Vectord> rDiag()const;

    /** \brief Diagonalization of a complex Matrix
    *
    * Uses Eigen to diagonalize the Matrix.
    * Phases eigenvectors s.t. first element is always positive.
    *
    * @param Transformation Matrix of eigenvectors
    * @param ev Vector of eigenvalues
    */
	void cDiag(Matrix<complex<double>>& Transformation, Vector<double>& ev) const;

    /** \brief Diagonalization of a complex Matrix
     *
     * Exactly like `Matrix::cDiag(Matrixcd&, Vectord&)` but returns copies instead of pass by reference
     *
     * @return Eigenvectors/value as a pair of <Matrixcd, Vectord>
     */
	pair<Matrix<complex<double>>, Vectord> cDiag()const;

	/** \brief Solve linear system of equations
	 *
	 * In the system of equations Ax=b, solve for x given A and b
	 *
	 * @param b RHS vector
	 * @return x Vector
	 */
	Vectord solveSLE(const Vectord& b);

	/** \brief Zero out the Matrix
	 *
	 */
	void zero();

	///@} End Math Operators

	//////////////////////////////////////////////////////////////////////
    /** \name I/O Operators
     * Group of functions for input and output
     */
    ///@{

	void print(ostream& os = cout) const;

	void write(const string& filename) const;

	virtual void write(ostream& os) const;

	void read(const string& filename);

	virtual void read(istream& os);

    ///@} End I/O Operators

    //////////////////////////////////////////////////////////////////////
    /** \name Getters/Setters */
    ///@{

	Vector<T> row(size_t i) const;
	Vector<T> col(size_t i) const;
	Vector<T> diag() const;

	size_t dim1() const { return dim1_; }

	size_t dim2() const { return dim2_; }

	T* ptr() const { return coeffs_; }

    ///@} End Getters/Setters

protected:
	T *coeffs_;
	size_t dim1_;
	size_t dim2_;
	bool ownership_;
};

/** \brief General typedef for complex matrices
 * \ingroup Core
 */
typedef Matrix<complex<double>> Matrixcd;

/** \brief General typedef for real matrices
 * \ingroup Core
 */
typedef Matrix<double> Matrixd;

/** \brief General typedef for Matrix<T>, Vectord pairs
 * \ingroup Core
 */
template <typename T>
using SpectralDecomposition = pair<Matrix<T>, Vectord>;

/** \brief Specialization of SpectralDecomposition of Matrixcd
 * \ingroup Core
 */
typedef SpectralDecomposition<complex<double>> SpectralDecompositioncd;

/** \brief Specialization of SpectralDecomposition of Matrixd
 * \ingroup Core
 */
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
Matrix<T> merge(const Matrix<T>& A, const Matrix<T>& B,
		const Matrix<T>& AB);

template<typename T>
Matrix<complex<double>> multATB(const Matrix<complex<double>>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> addAB(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> substAB(const Matrix<T>& A, const Matrix<T>& B);

template<typename T, typename U>
Matrix<T> multscalar(const U sca, const Matrix<T>& B);

template<typename T>
Matrix<T> re(const Matrix<T>& A);


//////////////////////////////////////////////////////////////////////
/// operator overloadings
//////////////////////////////////////////////////////////////////////

template <typename T>
Vector<T> operator*(const Matrix<T>& A, const Vector<T>& v);

//////////////////////////////////////////////////////////////////////
/// Diagonalization framework
//////////////////////////////////////////////////////////////////////

/**
 * \brief Diagonalize a complex double matrix
 * @param A Matrix to be diagonalized
 * @return Decomposed matrix
 */
SpectralDecompositioncd diagonalize(const Matrix<complex<double>>& A);

/**
 * \brief Diagonalize a complex double martrix
 * @param S Decomposed matrix (call-by-reference)
 * @param A Matrix to be diagonalized
 */
void diagonalize(SpectralDecompositioncd& S,const Matrix<complex<double>>& A);

/**
 * \brief Diagonalize a double matrix
 * @param A Matrix to be diagonalized
 * @return Decomposed matrix
 */
SpectralDecompositiond diagonalize(const Matrix<double>& A);

/**
 * \brief Diagonalize a double martrix
 * @param S Decomposed matrix (call-by-reference)
 * @param A Matrix to be diagonalized
 */
void diagonalize(SpectralDecompositiond& S,const Matrix<double>& A);

/**
 * \brief Diagonalization routine for remaining types (Eigenvector always double; not generally applicable)
 * @tparam T
 * @param A
 * @return
 */
template<typename T>
pair<Matrix<T>, Vectord> diagonalize(const Matrix<T>& A);

/**
 * \brief Construct matrix from its decomposition
 * @param X Decomposed matrix
 * @return re-constructed matrix
 */
template <typename T>
Matrix<T> toMatrix(const SpectralDecomposition<T>& X);

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
 * @return identity-matrix
 */
template<typename T>
Matrix<T> identityMatrix(size_t dim);

Matrixcd identityMatrixcd(size_t dim);
Matrixd identityMatrixd(size_t dim);

template<typename T>
Matrix<T> unitarySimilarityTrafo(const Matrix<T>& A,
		const Matrix<T>& B);

template<typename T>
Matrix<T> euclideanDistance(const Matrix<T>& A);

template<typename T>
double residual(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> regularize(const Matrix<T>& A, double eps);

template<typename T>
ostream& operator<<(ostream& os, const Matrix<T>& A);
template<typename T>
istream& operator>>(istream& is, Matrix<T>& A);

/// @TODO: Return Q & R matrix and adjust name convention
Matrixcd qr(const Matrixcd& A);

typedef tuple<Matrixcd, Matrixcd, Vectord> SVDcd;
typedef tuple<Matrixd, Matrixd, Vectord> SVDd;

SVDcd svd(const Matrixcd& A);
Matrixcd toMatrix(const SVDcd& svd);

SVDd svd(const Matrixd& A);

Eigen::MatrixXd toEigen(Matrixd A);
Eigen::MatrixXcd toEigen(Matrixcd A);
Matrixd toQutree(const Eigen::MatrixXd& A);
Matrixcd toQutree(const Eigen::MatrixXcd& A);

template <typename T>
Matrix<T> subMatrix(const Matrix<T> A, size_t dim1, size_t dim2);

/**
 * @TODOs:
 * - move non-BLAS functions to matrix extension (Christoph)
 * - matrix inherits from vector (Christoph)
 * - general elementwise function for Vectors, matrices & tensors (Roman)
 * - range operators for Vector & Matrix class (Roman)
 * - scalar functions for Matrices (f(A)) sqrt(x), matrixfunction(A, sqrt); (Roman)
 * - DiagonalMatrix class (inherit from Vector): add and multiply member-functions
 * - Hermitisch & unitär-inheritance
 * - extend Unit tests
 **/

