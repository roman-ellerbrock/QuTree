#pragma once
#include "stdafx.h"

/**
 * \class Vector
 * \ingroup QD-lib
 * \brief This class is a Vector represented in a basis.
 *
 * The class is a simple version of a represented Vector.
 * The class is suited for Quantum Dynamics simulations; it is
 * not a generell purpose library class.
 *
 * */

template <typename T>
class Vector
{
public:
	//////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////
	// Standard Constructor
	Vector():Vector(1) {}

	// Constructor
	Vector(size_t dim_);

	// Copy constructor
	Vector(const Vector& old);

	// Move constructor
	Vector(Vector&& old);

	// Copy Assignment Operator
	Vector& operator=(const Vector& other);

	// Move Assignment Operator
	Vector& operator=(Vector&& other);

	// Destructor
	~Vector();

	// I/O
	void Write(ostream& os);

	void print(ostream& os = cout)const;

	// Operators
	inline T operator()(size_t i)const
	{
		assert(i < dim && i >= 0);
		return coeffs[i];
	}

	inline T& operator()(size_t i)
	{
		assert(i < dim && i >= 0);
		return coeffs[i];
	}

	Vector operator+(const Vector b);
	Vector operator-(const Vector b);
	T operator*(const Vector b);
	Vector& operator*=(T coeff);

	// Math Operators
	double Norm()const;

	void zero();

	// Setter & Getter
	inline size_t Dim()const { return dim; }

protected:
	T* coeffs;
	size_t dim;

};

typedef Vector<double> Vectord;
typedef Vector<complex<double>> Vectorcd;
