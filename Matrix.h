#pragma once

#include "stdafx.h"
#include <Eigen/Dense>
#include "Vector.h"

/**
 * \class Matrix
 * \ingroup QD-lib
 * \brief This class represents a Matrix.
 *
 * This is a simple Matrix-class. The functions that
 * are implemented are specially suited for QuantumDynamics
 * simulations.
 *
 * */


// @TODO: Check constructors and =-operators

template<typename T>
class Matrix {
public:
	//////////////////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////////////////

	Matrix();

	Matrix(size_t dim1_, size_t dim2_);

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
	// Fundamental Management
	//////////////////////////////////////////////////////////////////////

    T& operator()(const size_t i, const size_t j)const;

	T& operator()(const size_t i, const size_t j);

	inline T operator[](const size_t idx) const {
		return coeffs[idx];
	}

	inline T& operator[](const size_t idx) {
		return coeffs[idx];
	}

	//////////////////////////////////////////////////////////////////////
	// Fundamental Math operators
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

	//////////////////////////////////////////////////////////////////////
	// More Math operators
	//////////////////////////////////////////////////////////////////////

	// Calculate the Trace
	T Trace() const;

	// Transpose a matrix
	Matrix Transpose();

	// Inverts a complex matrix
	Matrix<complex<double> > cInv() const;

	// Diagonalization
	void rDiag(Matrix<double>& Transformation, Vector<double>& ev);

	void cDiag(Matrix<complex<double>>& Transformation, Vector<double>& ev) const;

	// Solve System of linear equations for Matrix<double>
	Vectord SolveSLE(const Vectord& b_);

	void Zero();

	//////////////////////////////////////////////////////////////////////
	// I/O
	//////////////////////////////////////////////////////////////////////

	void print(ostream& os = cout) const;

	void Write(ostream& os) const;

	void Read(istream& os);


	//////////////////////////////////////////////////////////////////////
	// Getter & Setter
	//////////////////////////////////////////////////////////////////////

	size_t Dim1() const { return dim1; }

	size_t Dim2() const { return dim2; }

	size_t Size() const { return size; }

	T *Coeffs() const { return coeffs; }

protected:
	double conjugate(const double d) const {
		return d;
	}

	complex<double> conjugate(const complex<double> c) const {
		return conj(c);
	}

	T *coeffs;
	size_t dim1;
	size_t dim2;
	size_t size;
};

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

template<typename T>
Matrix<T> UnitarySimilarityTrafo(const Matrix<T>& A,
		const Matrix<T>& B);

template<typename T>
Matrix<T> Regularize(const Matrix<T>& A, double eps);

typedef Matrix<complex<double>> Matrixcd;
typedef Matrix<double> Matrixd;

