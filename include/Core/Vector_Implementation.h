#pragma once
#include "Vector.h"

//////////////////////////////////////////////////////////////////////
/// Constructors, Destructors, Asignment operators
/// ~ Rule of Five ~
//////////////////////////////////////////////////////////////////////
// Constructor
template<typename T>
Vector<T>::Vector(size_t dim)
	:dim_(dim), coeffs_(new T[dim]) {
	assert(dim > 0);
	zero();
}

template<typename T>
Vector<T>::Vector(const string& filename)
	:Vector() {
	read(filename);
}

template<typename T>
Vector<T>::Vector(istream& is)
	:Vector() {
	read(is);
}

// Copy constructor
template<typename T>
Vector<T>::Vector(const Vector& old)
	:Vector(old.dim_) {
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] = old.coeffs_[i];
}

// Move constructor
template<typename T>
Vector<T>::Vector(Vector&& old) noexcept
	:dim_(old.dim_), coeffs_(old.coeffs_) {
	old.coeffs_ = nullptr;
}

// Copy Assignment Operator
template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
	Vector tmp(other);
	*this = move(tmp);
	return *this;
}

// Move Assignment Operator
template<typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& other) noexcept {
	delete[] coeffs_;
	coeffs_ = other.coeffs_;
	other.coeffs_ = nullptr;
	dim_ = other.dim_;
	return *this;
}

// Destructos
template<typename T>
Vector<T>::~Vector() {
	delete[] coeffs_;
}

//////////////////////////////////////////////////////////////////////
// File I/O
//////////////////////////////////////////////////////////////////////
template<typename T>
void Vector<T>::write(ostream& os) const {
	// Verification
	os.write("VECT", 4);
	os.write((char *) &dim_, sizeof(dim_));

	int32_t size = sizeof(T);
	os.write((char *) &size, sizeof(size));
	for (size_t i = 0; i < dim_; i++) {
		T Coeff_now = operator()(i);
		os.write((char *) &Coeff_now, size);
	}
	os.flush();
}

template<typename T>
void Vector<T>::read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("VECT");
	assert(s_key == s_check);

	// Read dimensions
	is.read((char *) &dim_, sizeof(dim_));
	(*this) = Vector<T>(dim_);

	// Read the size of type
	int32_t size;
	is.read((char *) &size, sizeof(size));
	assert(size == sizeof(T));

	// Read the coefficients
	for (size_t i = 0; i < dim_; i++) {
		T Coeff_now;
		is.read((char *) &Coeff_now, size);
		operator[](i) = Coeff_now;
	}
}

template<typename T>
void Vector<T>::write(const string& filename) const {
	ofstream os(filename);
	write(os);
}

template<typename T>
void Vector<T>::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

template<typename T>
void Vector<T>::print(ostream& os) const {
	for (size_t i = 0; i < dim_; i++)
		os << coeffs_[i] << " ";
	os << endl;
}

//////////////////////////////////////////////////////////////////////
// Arithmetic
//////////////////////////////////////////////////////////////////////
template<typename T>
Vector<T> Vector<T>::operator+=(Vector b) {
	assert(b.dim() == dim());
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] += b[i];
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator-=(Vector b) {
	assert(b.dim() == dim());
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] -= b[i];
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T> b) const {
	Vector res(*this);
	res += b;
	return res;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T> b) const {
	Vector res(*this);
	res -= b;
	return res;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(T coeff) {
	for (size_t i = 0; i < dim_; i++) {
		coeffs_[i] *= coeff;
	}
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator/=(T coeff) {
	for (size_t i = 0; i < dim_; i++) {
		coeffs_[i] /= coeff;
	}
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator*(T coeff) const {
	Vector<T> C(*this);
	C *= coeff;
	return C;
}

template<typename T>
Vector<T> Vector<T>::operator/(T coeff) const {
	Vector<T> C(*this);
	C /= coeff;
	return C;
}

template<typename T>
T Vector<T>::operator*(const Vector<T> b) const {
	assert(b.dim() == dim_);
	T res = 0;
	for (size_t i = 0; i < dim_; i++)
		res += conj(operator()(i)) * b(i);
	return res;
}

template<typename T>
double Vector<T>::norm() const {
	double norm = 0;
	for (size_t i = 0; i < dim_; i++) {
		norm += pow(abs(operator()(i)), 2);
	}
	norm = sqrt(norm);
	return norm;
}

template<typename T>
void Vector<T>::zero() {
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] = 0;
}

template<typename T>
T Sum(Vector<T>& a) {
	T sum = 0.;
	for (size_t i = 0; i < a.dim(); ++i) {
		sum += a(i);
	}
	return sum;
}

template <typename T>
void normalize(Vector<T>& a) {
	double norm = a.norm();
	a /= (norm);
}

template<typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b) {
	double distance = 0;
	for (size_t i = 0; i < a.dim(); ++i) {
		T d = a(i) - b(i);
		distance += pow(abs(d), 2);
	}
	return sqrt(distance);
}

template<typename T>
double residual(const Vector<T>& A, const Vector<T>& B) {
	return euclidean_distance(A, B);
}

template<typename T>
Vector<T> reverse(const Vector<T>& A) {
	Vector<T> B(A.dim());
	for (size_t i = 0; i < A.dim(); ++i) {
		B(A.dim() - 1 - i) = A(i);
	}
	return B;
}

template<typename T>
T sum(const Vector<T>& A) {
	T s = 0.;
	for (size_t i = 0; i < A.dim(); ++i) {
		s+= A(i);
	}
	return s;
}

template<typename T>
Vector<T> regularize(Vector<T> A, double eps) {
	for (size_t i = 0; i < A.dim(); ++i) {
		A(i) += eps * exp(-A(i) / eps);
	}
	return A;
}

template<typename T>
Vector<T> inverse(Vector<T> A, double eps) {
	A = regularize(A, eps);
	for (size_t i = 0; i < A.dim(); ++i) {
		A(i) = 1. / A(i);
	}
	return A;
}

Vectord toQutree(const Eigen::VectorXd& v) {
	Vectord vqutree(v.rows());
	for (size_t i = 0; i < vqutree.dim(); ++i) {
		vqutree(i) = v(i);
	}
	return vqutree;
}

Vectorcd toQutree(const Eigen::VectorXcd& v) {
	Vectorcd vqutree(v.rows());
	for (size_t i = 0; i < vqutree.dim(); ++i) {
		vqutree(i) = v(i);
	}
	return vqutree;
}

template<typename T>
void elementwise(Vector<T>& res, const Vector<T>& A, const function<T(T)>& f) {
	assert(A.dim() == res.dim());
	for (size_t i = 0; i < A.Dim1()*A.Dim2(); ++i) {
		res[i] = f(A[i]);
	}
}

template <typename T>
Vector<T> elementwise(const Vector<T>& A, const function<T(T)>& f) {
	Vector<T> res(A.dim());
	elementwise(res, A, f);
	return res;
}
