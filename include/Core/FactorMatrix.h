#pragma once
#include "Matrix.h"
#include "Tensor.h"


template<typename T>
class FactorMatrix: public Matrix<T> {
/**
 * \class FactorMatrix
 * \ingroup Core
 * \brief This class represents a factor matrix.
 *
 * Factor Matrices can be multiplied with tensors and result from
 * performing tensor-hole products.
 */
public:
	FactorMatrix() : mode_(0) {}

	FactorMatrix(const Matrix<T>& A, size_t k)
		: Matrix<T>(A), mode_(k) {
	}

	explicit FactorMatrix(size_t dim_a, size_t dim_b, size_t mode)
		: Matrix<T>(dim_a, dim_b),  mode_(mode) {}

	FactorMatrix(size_t dim_, size_t mode)
		: Matrix<T>(dim_, dim_), mode_(mode) {}

	FactorMatrix(const TensorDim& tdim, size_t mode)
		: Matrix<T>(tdim.Active(mode), tdim.Active(mode)), mode_(mode) {
	}

	explicit FactorMatrix(istream& is) {
		Read(is);
	}

	explicit FactorMatrix(const string& filename) {
		ifstream is(filename);
		Read(is);
	}

	/////////////////////////////////////////////////////
	// Rule of five
	/////////////////////////////////////////////////////

	// Copy constructor
	FactorMatrix(const FactorMatrix& A)
		: Matrix<T>(A), mode_(A.Mode()) {
	}

	// Move constructor
	FactorMatrix(FactorMatrix&& A) noexcept
		: Matrix<T>(A), mode_(A.Mode()) {
	}

	// copy assignment
	FactorMatrix& operator=(const FactorMatrix& other) {
		Matrix<T>::operator=(other);
		mode_ = other.mode_;
		return *this;
	}

	// move assignment
	FactorMatrix& operator=(FactorMatrix&& other) noexcept {
		Matrix<T>::operator=(other);
		mode_ = other.mode_;
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

	size_t Mode() const { return mode_; }
	size_t& Mode() { return mode_; }

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

	void Write(ostream& os) const override {
		os.write("FAMA", 4);
		os.write((char *) &mode_, sizeof(mode_));
		Matrix<T>::Write(os);
	}

	void Write(const string& filename) const {
		ofstream os(filename);
		Write(os);
	}

	void Read(istream& is) override {
		char check[5];
		is.read(check, 4);
		string s_check(check, 4);
		string s_key("FAMA");
		assert(s_key == s_check);
		is.read((char *) &mode_, sizeof(mode_));
		Matrix<T>::Read(is);
	}

	void Read(const string& filename) {
		ifstream is(filename);
		Read(is);
	}

protected:
	size_t mode_;
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
FactorMatrix<T> HoleProduct(const Tensor<T>& A, const Tensor<T>& B, size_t mode) {
	const TensorDim& tdim_a = A.Dim();
	const TensorDim& tdim_b = B.Dim();
	size_t act_a = tdim_a.Active(mode);
	size_t act_b = tdim_b.Active(mode);
	FactorMatrix<T> S(act_a, act_b, mode);
	if (mode == tdim_a.GetOrder()) {
		cout << "active : " << act_a << endl;
		cout << "active : " << act_b << endl;
	}
	HoleProduct(S, A, B, mode);
	return S;
}

template<typename T>
void HoleProduct(FactorMatrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t mode) {
	const TensorDim& tdim_a = A.Dim();
	size_t bef = tdim_a.Before(mode);
	size_t act_a = tdim_a.Active(mode);
	size_t aft = tdim_a.TotAfter(mode);
	const TensorDim& tdim_b = B.Dim();
	size_t act_b = tdim_b.Active(mode);
	assert(S.Dim1() == act_a);
	assert(S.Dim2() == act_b);
	TensorHoleProduct(S, A, B, bef, act_a, act_b, aft);
}

template<typename T>
void DotProduct(FactorMatrix<T>& S, const Tensor<T>& A, const Tensor<T>& B) {
	const TensorDim& tdim = A.Dim();
	size_t mode = tdim.GetOrder();
	size_t bef = tdim.Before(mode);
	size_t act = tdim.Active(mode);
	size_t aft = tdim.TotAfter(mode);
	TensorHoleProduct(S, A, B, bef, act, act, aft);
}

template<typename T>
FactorMatrix<T>
SPOUnitarySimilarityTrafo(const FactorMatrix<T>& A, const FactorMatrix<T>& B) {
	assert(A.Mode() == B.Mode());
	Matrix<T> resultmat = UnitarySimilarityTrafo(A, B);
	return FactorMatrix<T>(resultmat, A.Mode());
}

template<typename T>
ostream& operator<<(ostream& os, const FactorMatrix<T>& A) {
	A.Write(os);
	return os;
}

template<typename T>
istream& operator>>(istream& is, FactorMatrix<T>& A) {
	A.Read(is);
	return is;
}

typedef FactorMatrix<double> FactorMatrixd;

typedef FactorMatrix<complex<double>> FactorMatrixcd;
