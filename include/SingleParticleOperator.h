#pragma once
#include "Matrix.h"
#include "Tensor.h"


template<typename T>
class SingleParticleOperator:
	public Matrix<T> {
public:
	SingleParticleOperator() = default;

	SingleParticleOperator(const Matrix<T>& A, size_t k)
		: Matrix<T>(A), mode(k) {
	}

	SingleParticleOperator(size_t dim_, size_t mode_)
		: Matrix<T>(dim_, dim_), mode(mode_) {}

	SingleParticleOperator(const TensorDim& tdim, size_t mode_)
		: Matrix<T>(tdim.Active(mode_), tdim.Active(mode_)), mode(mode_) {}

	/////////////////////////////////////////////////////
	// Rule of five
	/////////////////////////////////////////////////////

	// Copy constructor
	SingleParticleOperator(const SingleParticleOperator& A)
		: Matrix<T>(A), mode(A.Mode()) {
	}

	// Move constructor
	SingleParticleOperator(SingleParticleOperator&& A) noexcept
		: Matrix<T>(A), mode(A.Mode()) {
	}

	// copy assignment
	SingleParticleOperator& operator=(const SingleParticleOperator& other) {
		Matrix<T>::operator=(other);
		mode = other.mode;
		return *this;
	}

	// move assignment
	SingleParticleOperator& operator=(SingleParticleOperator&& other) noexcept {
		Matrix<T>::operator=(other);
		mode = other.mode;
		return *this;
	}

	~SingleParticleOperator() = default;

	SingleParticleOperator& operator+=(const SingleParticleOperator& other) {
		assert(other.Mode() == this->Mode());
		const Matrix<T>& other_mat(other);
		return (SingleParticleOperator&) Matrix<T>::operator+=(other_mat);
	}

	/////////////////////////////////////////////////////
	// Getter & Setter
	/////////////////////////////////////////////////////

	size_t Dim() const { return Matrix<T>::Dim1(); }

	size_t Mode() const { return mode; }

	// TensorC = A * TensorB
	template<typename U>
	friend Tensor<U> operator*(const SingleParticleOperator<T>& A,
		const Tensor<U>& B) {
		return multAB(A, B);
	}

	// C = A * B
	friend SingleParticleOperator<T> operator*(const SingleParticleOperator<T>& A,
		const SingleParticleOperator<T>& B) {
		Matrix<T> C = Matrix<T>(A) * Matrix<T>(B);
		return SingleParticleOperator<T>(C, B.Mode());
	}

	// C = A + B
	friend SingleParticleOperator<T> operator+(const SingleParticleOperator<T>& A,
		const SingleParticleOperator<T>& B) {
		Matrix<T> C = Matrix<T>(A) + Matrix<T>(B);
		return SingleParticleOperator<T>(C, B.Mode());
	}

	// B = scalar * A
	template<typename U>
	friend SingleParticleOperator<T> operator*(U scalar,
		const SingleParticleOperator<T>& A) {
		Matrix<T> B = scalar * Matrix<T>(A);
		return SingleParticleOperator<T>(B, A.Mode());
	}

protected:
	size_t mode;
};

template<typename T, typename U>
void multAB(Tensor<T>& C, const SingleParticleOperator<U>& A, const Tensor<T>& B, bool zero = true) {
	int mode = A.Mode();
	return multAB(C, A, B, mode, zero);
}

template<typename T, typename U>
Tensor<T> multAB(const SingleParticleOperator<U>& A, const Tensor<T>& B) {
	int mode = A.Mode();
	return multAB(A, B, mode);
}

template<typename T, typename U>
Tensor<T> multATB(const SingleParticleOperator<U>& A, const Tensor<T>& B) {
	int mode = A.Mode();
	return multATB(A, B, mode);
}

template<typename T>
SingleParticleOperator<T>
SPOUnitarySimilarityTrafo(const SingleParticleOperator<T>& A, const SingleParticleOperator<T>& B) {
	assert(A.Mode() == B.Mode());
	Matrix<T> resultmat = UnitarySimilarityTrafo(A, B);
	return SingleParticleOperator<T>(resultmat, A.Mode());
}

typedef SingleParticleOperator<double> SPOd;

typedef SingleParticleOperator<complex<double>> SPOcd;
