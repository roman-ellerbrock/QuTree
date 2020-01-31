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
		: Vector(1) {}

	// Constructor
	explicit Vector(size_t dim_);

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
		assert(i < dim);
		return coeffs[i];
	}

	inline T& operator()(size_t i) {
		assert(i < dim);
		return coeffs[i];
	}

	inline T operator[](size_t i) const {
		return coeffs[i];
	}

	inline T& operator[](size_t i) {
		return coeffs[i];
	}

	Vector operator+(Vector b);
	Vector operator-(Vector b);
	T operator*(Vector b);
	Vector& operator*=(T coeff);
	Vector& operator/=(T coeff);

	// Math Operators
	double Norm() const;

	void zero();

	// Setter & Getter
	inline size_t Dim() const { return dim; }

protected:
	T *coeffs;
	size_t dim;
};

template<typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b);

template<typename T>
double Residual(const Vector<T>& A, const Vector<T>& B);

typedef Vector<double> Vectord;

typedef Vector<complex<double>> Vectorcd;
