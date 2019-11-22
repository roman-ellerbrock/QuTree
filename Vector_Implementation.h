#pragma once
#include "Vector.h"

// Constructor
template <typename T>
Vector<T>::Vector(size_t dim_)
	:dim(dim_),	coeffs(new T[dim_])
{
	assert(dim > 0);
	zero();
}

// Copy constructor
template <typename T>
Vector<T>::Vector(const Vector& old)
	:Vector(old.dim)
{
	for (size_t i = 0; i < dim; i++)
		coeffs[i] = old.coeffs[i];
}

// Move constructor
template <typename T>
Vector<T>::Vector(Vector&& old)
	:dim(old.dim), coeffs(old.coeffs)
{
	old.coeffs = nullptr;
}

// Copy Assignment Operator
template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other)
{
	Vector tmp(other);
	*this = move(tmp);
	return *this;
}

// Move Assignment Operator
template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& other)
{
	delete[] coeffs;
	coeffs = other.coeffs;
	other.coeffs = nullptr;
	dim = other.dim;
	return *this;
}

// Destructos
template <typename T>
Vector<T>::~Vector()
{
	delete[] coeffs;
}

// Math Operators
template <typename T>
void Vector<T>::Write(ostream& os)
{
	for (size_t i = 0; i < dim; i++)
	{
		os << operator()(i) << endl;
	}
}

template <typename T>
void Vector<T>::print(ostream& os)const
{
	for (size_t i = 0; i < dim; i++)
		os << coeffs[i] << " ";
	os << endl;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T> b)
{
	assert(b.Dim() == dim);
	Vector res(dim);
	for (size_t i = 0; i < dim ; i++)
		res(i) = operator()(i) + b(i);
	return res;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T> b)
{
	assert(b.Dim() == dim);
	Vector res(dim);
	for (size_t i = 0; i < dim ; i++)
		res(i) = operator()(i) - b(i);
	return res;
}

template <typename T>
T Vector<T>::operator*(const Vector<T> b)
{
	assert(b.Dim() == dim);
	T res = 0;
	for (size_t i = 0; i < dim ; i++)
		res += operator()(i) * b(i);
	return res;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(T coeff) {
	for (size_t i = 0; i < dim; i++){
		coeffs[i] *= coeff;
	}
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator/=(T coeff) {
	for (size_t i = 0; i < dim; i++){
		coeffs[i] /= coeff;
	}
	return *this;
}

// Math Operators
template <typename T>
double Vector<T>::Norm()const
{
	double norm = 0;

	for (size_t i = 0; i < dim; i++)
	{
		norm += pow(abs(operator()(i)), 2);
	}

	norm = sqrt(norm);
	
	return norm;
}

template <typename T>
void Vector<T>::zero()
{
	for (size_t i = 0; i < dim; i++)
		coeffs[i] = 0;
}

template <typename T>
double euclidean_distance(const Vector<T>& a, const Vector<T>& b) {
	double distance = 0;
	for (size_t i = 0; i < a.Dim(); ++i) {
		distance += pow(a(i) - b(i), 2);
	}
	return sqrt(distance);
}

