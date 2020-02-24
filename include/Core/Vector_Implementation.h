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
	Zero();
}

template<typename T>
Vector<T>::Vector(const string& filename)
	:Vector() {
	Read(filename);
}

template<typename T>
Vector<T>::Vector(istream& is)
	:Vector() {
	Read(is);
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
void Vector<T>::Write(ostream& os) const {
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
void Vector<T>::Read(istream& is) {
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
void Vector<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template<typename T>
void Vector<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
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
Vector<T> Vector<T>::operator+(const Vector<T> b) {
	assert(b.Dim() == dim_);
	Vector res(dim_);
	for (size_t i = 0; i < dim_; i++)
		res(i) = operator()(i) + b(i);
	return res;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T> b) {
	assert(b.Dim() == dim_);
	Vector res(dim_);
	for (size_t i = 0; i < dim_; i++)
		res(i) = operator()(i) - b(i);
	return res;
}

template<typename T>
Vector<T> Vector<T>::operator+=(Vector b) {
	assert(b.Dim() == Dim());
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] += b[i];
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator-=(Vector b) {
	assert(b.Dim() == Dim());
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] -= b[i];
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator*(T coeff) {
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] *= coeff;
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator/(T coeff) {
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] /= coeff;
	return *this;
}

template<typename T>
T Vector<T>::operator*(const Vector<T> b) {
	assert(b.Dim() == dim_);
	T res = 0;
	for (size_t i = 0; i < dim_; i++)
		res += operator()(i) * b(i);
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
double Vector<T>::Norm() const {
	double norm = 0;
	for (size_t i = 0; i < dim_; i++) {
		norm += pow(abs(operator()(i)), 2);
	}
	norm = sqrt(norm);
	return norm;
}

template<typename T>
void Vector<T>::Zero() {
	for (size_t i = 0; i < dim_; i++)
		coeffs_[i] = 0;
}

template<typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b) {
	double distance = 0;
	for (size_t i = 0; i < a.Dim(); ++i) {
		T d = a(i) - b(i);
		distance += pow(abs(d), 2);
	}
	return sqrt(distance);
}

template<typename T>
double Residual(const Vector<T>& A, const Vector<T>& B) {
	return euclidean_distance(A, B);
}

template<typename T>
Vector<T> Inverse(Vector<T> A, double eps) {
	for (size_t i = 0; i < A.Dim(); ++i) {
		A(i) = 1. / (A(i) + eps * exp(-A(i) / (A(i) + eps)));
	}
	return A;
}

