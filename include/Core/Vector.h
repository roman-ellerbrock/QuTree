#pragma once
#include "stdafx.h"
#include "Eigen/Dense"

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
	void write(const string& filename) const;
	void write(ostream& os) const;

	void read(const string& filename);
	void read(istream& is);

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

	friend Vector operator*(T a, const Vector<T>& A) {
		Vector B(A);
		B *= a;
		return B;
	}

	// Math Operators
	double norm() const;

	void zero();

	// Setter & Getter
	inline size_t dim() const { return dim_; }

protected:
	T *coeffs_;
	size_t dim_;
};

template <typename T>
void normalize(Vector<T>& a);

template<typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b);

template<typename T>
double residual(const Vector<T>& A, const Vector<T>& B);

template<typename T>
Vector<T> reverse(const Vector<T>& A);

template<typename T>
T sum(const Vector<T>& A);

template<typename T>
Vector<T> regularize(Vector<T> A, double eps);

template<typename T>
Vector<T> inverse(Vector<T> A, double eps = 1e-7);

typedef Vector<double> Vectord;

typedef Vector<complex<double>> Vectorcd;

Vectord toQutree(const Eigen::VectorXd& v);

Vectorcd toQutree(const Eigen::VectorXcd& v);

template <typename T>
void elementwise(Vector<T>& res, const Vector<T>& A, const function<T(T)>& f);

template <typename T>
Vector<T> elementwise(const Vector<T>& A, const function<T(T)>& f);

