#pragma once
#include "Matrix.h"
#include "Tensor.h"


template<typename T>
class FactorMatrix: public Matrix<T> {
public:
	FactorMatrix() : mode(0) {}

	FactorMatrix(const Matrix<T>& A, size_t k)
		: Matrix<T>(A), mode(k) {
	}

	FactorMatrix(size_t dim_, size_t mode_)
		: Matrix<T>(dim_, dim_), mode(mode_) {}

	FactorMatrix(const TensorDim& tdim, size_t mode_)
		: Matrix<T>(tdim.Active(mode_), tdim.Active(mode_)), mode(mode_) {}

	/////////////////////////////////////////////////////
	// Rule of five
	/////////////////////////////////////////////////////

	// Copy constructor
	FactorMatrix(const FactorMatrix& A)
		: Matrix<T>(A), mode(A.Mode()) {
	}

	// Move constructor
	FactorMatrix(FactorMatrix&& A) noexcept
		: Matrix<T>(A), mode(A.Mode()) {
	}

	// copy assignment
	FactorMatrix& operator=(const FactorMatrix& other) {
		Matrix<T>::operator=(other);
		mode = other.mode;
		return *this;
	}

	// move assignment
	FactorMatrix& operator=(FactorMatrix&& other) noexcept {
		Matrix<T>::operator=(other);
		mode = other.mode;
		return *this;
	}

	~FactorMatrix() = default;

	FactorMatrix& operator+=(const FactorMatrix& other) {
		assert(other.Mode() == this->Mode());
		const Matrix<T>& other_mat(other);
		return (FactorMatrix&) Matrix<T>::operator+=(other_mat);
	}

	/////////////////////////////////////////////////////
	// Getter & Setter
	/////////////////////////////////////////////////////

	size_t Dim() const { return Matrix<T>::Dim1(); }

	size_t Mode() const { return mode; }

	// TensorC = A * TensorB
	template<typename U>
	friend Tensor<U> operator*(const FactorMatrix<T>& A,
		const Tensor<U>& B) {
		return multAB(A, B);
	}

	// C = A * B
	friend FactorMatrix<T> operator*(const FactorMatrix<T>& A,
		const FactorMatrix<T>& B) {
		Matrix<T> C = Matrix<T>(A) * Matrix<T>(B);
		return FactorMatrix<T>(C, B.Mode());
	}

	// C = A + B
	friend FactorMatrix<T> operator+(const FactorMatrix<T>& A,
		const FactorMatrix<T>& B) {
		Matrix<T> C = Matrix<T>(A) + Matrix<T>(B);
		return FactorMatrix<T>(C, B.Mode());
	}

	// B = scalar * A
	template<typename U>
	friend FactorMatrix<T> operator*(U scalar,
		const FactorMatrix<T>& A) {
		Matrix<T> B = scalar * Matrix<T>(A);
		return FactorMatrix<T>(B, A.Mode());
	}

protected:
	size_t mode;
};

template<typename T, typename U>
void multAB(Tensor<T>& C, const FactorMatrix<U>& A, const Tensor<T>& B, bool zero = true) {
	int mode = A.Mode();
	return multAB(C, A, B, mode, zero);
}

template<typename T, typename U>
Tensor<T> multAB(const FactorMatrix<U>& A, const Tensor<T>& B) {
	int mode = A.Mode();
	return multAB(A, B, mode);
}

template<typename T, typename U>
Tensor<T> multATB(const FactorMatrix<U>& A, const Tensor<T>& B) {
	int mode = A.Mode();
	return multATB(A, B, mode);
}

template<typename T>
FactorMatrix<T>
SPOUnitarySimilarityTrafo(const FactorMatrix<T>& A, const FactorMatrix<T>& B) {
	assert(A.Mode() == B.Mode());
	Matrix<T> resultmat = UnitarySimilarityTrafo(A, B);
	return FactorMatrix<T>(resultmat, A.Mode());
}

typedef FactorMatrix<double> FactorMatrixd;

typedef FactorMatrix<complex<double>> FactorMatrixcd;
