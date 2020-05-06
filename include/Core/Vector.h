#pragma once
#include "stdafx.h"


template<typename T>
class Vector {
/**
 * \class Vector
 * \ingroup Core
 * \brief This class is a Vector represented in a basis.
 *
 * The class is a simple version of a represented Vector.
 * The class is suited for Quantum Dynamics simulations; it is
 * not a generell purpose library class.
 *
 * */
public:
	//////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////
	// Standard Constructor
	Vector()
		: Vector(1)  {}

	// Constructor
	explicit Vector(size_t dim);

	explicit Vector(const string& filename);

	explicit Vector(istream& is);

	// Copy constructor
	Vector(const Vector& old);

	// Move constructor
	Vector(Vector&& old) noexcept;

	// Copy Assignment Operator
	Vector& operator=(const Vector& other);

	// Move Assignment Operator
	Vector& operator=(Vector&& other) noexcept;

	// Destructor
	~Vector();

	// I/O
	void Write(const string& filename) const;
	void Write(ostream& os) const;

	void Read(const string& filename);
	void Read(istream& is);

	void print(ostream& os = cout) const;

	// Operators
	inline T operator()(size_t i) const {
		assert(i < dim_);
		return coeffs_[i];
	}

	inline T& operator()(size_t i) {
		assert(i < dim_);
		return coeffs_[i];
	}

	inline T operator[](size_t i) const {
		return coeffs_[i];
	}

	inline T& operator[](size_t i) {
		return coeffs_[i];
	}

	Vector operator+=(Vector b);
	Vector operator-=(Vector b);
	Vector operator+(Vector b)const;
	Vector operator-(Vector b)const;
	T operator*(Vector b)const;
	Vector& operator*=(T coeff);
	Vector& operator/=(T coeff);
	Vector operator*(T coeff)const;
	Vector operator/(T coeff)const;

	// Math Operators
	double Norm() const;

	void Zero();

	// Setter & Getter
	inline size_t Dim() const { return dim_; }

protected:
	T *coeffs_;
	size_t dim_;
};

template<typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b);

template<typename T>
double Residual(const Vector<T>& A, const Vector<T>& B);

template<typename T>
Vector<T> reverse(const Vector<T>& A);

template<typename T>
T sum(const Vector<T>& A);

template<typename T>
Vector<T> Regularize(Vector<T> A, double eps);

template<typename T>
Vector<T> Inverse(Vector<T> A, double eps = 1e-7);

typedef Vector<double> Vectord;

typedef Vector<complex<double>> Vectorcd;
