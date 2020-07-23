#pragma once
#include "TensorShape.h"
#include "Core/Matrix.h"

/**
 * \defgroup Core
 * \brief This group includes the basic datastructures in QuTree.
 *
 * These datastructures include the Vector, Matrix, and Tensor classes.
 */

template <typename T>
class Tensor
/**
 * \class Tensor
 * \ingroup Core
 * \brief This class represents a set of mathematical Tensors of n-th order.
 *
 * The class allows to handle arithmentic operations on Tensors like
 * matrix-products, etc. The set of tensors can be interpreted simultaneously
 * as a set of vectors and thus operations are available like orthogonalizations,
 * Dot-products, etc. Superindex mappings are used frequently throughout the
 * class functions.
 *
 * Usage:
 * TensorDim dim_({2, 3, 4}, 1);
 * Tensorcd A(dim_);
 *
 * */
{
public:
	//////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////

	// Standard Constructor
	Tensor() : coeffs_(new T[1]) {}

	Tensor(const initializer_list<size_t>& dim, bool InitZero = true);

	// Constructor with TensorDim
	explicit Tensor(const TensorShape& dim, bool InitZero = true);

	// Construct from external memory
	explicit Tensor(const TensorShape& dim, T* ptr, bool ownership = true, bool InitZero = true);

	explicit Tensor(istream& is);

	explicit Tensor(const string& filename);

	// Construct Tensor from Matrix
	explicit Tensor(const Matrix<T>& mat);

	// Copy constructor
	Tensor(const Tensor& old);

	// Copy-Multiply constructor
	Tensor(const Tensor& old, T factor);

	// Move constructor
	Tensor(Tensor&& old)noexcept;

	// Copy Assignment Operator
	Tensor& operator=(const Tensor& old);

	// Move Assignment Operator
	Tensor& operator=(Tensor&& old)noexcept;

	// Destructor
	~Tensor();

	//////////////////////////////////////////////////////////
	// File handling
	//////////////////////////////////////////////////////////
	void print(ostream& os = cout)const;

	void Write(ostream& os)const;

	void Write(const string& filename)const;

	void Read(istream& is);

	void Read(const string& filename);

	//////////////////////////////////////////////////////////
	// Bracket Operators
	//////////////////////////////////////////////////////////
	T& operator()(size_t i)const;

	T& operator()(size_t i);

	const T& operator()(size_t i, size_t n)const;

	T& operator()(size_t i, size_t n);

	T& operator()(size_t bef, size_t j, size_t aft, size_t leaf);

	const T& operator()(size_t bef, size_t j, size_t aft, size_t leaf)const;

	T& operator()(size_t i, size_t j, size_t k, size_t f, size_t n);

	const T& operator()(size_t i, size_t j, size_t k, size_t f, size_t n)const;

	T& operator()(const vector<size_t>& dims);

	const T& operator()(const vector<size_t>& dims) const;

	// double hole operator
	T& operator()(size_t bef, size_t i, size_t mid, size_t j, size_t beh,
		size_t mode1, size_t mode2, size_t n)const;

	inline const T& operator[](const size_t idx)const
	{
		// Fast bracket operator
		return coeffs_[idx];
	}

	inline T& operator[](const size_t idx)
	{
		// Fast bracket operator
		return coeffs_[idx];
	}

	//////////////////////////////////////////////////////////
	// Math Operators
	//////////////////////////////////////////////////////////
	friend Tensor<T> operator+(const Tensor<T>& A, const Tensor<T>& B)
	{
		assert(A.shape().GetDimTot() == B.shape().GetDimTot());
		Tensor C(A.shape());
		for (int i = 0; i < A.shape().GetDimTot(); i++)
		{
			C(i) = A(i) + B(i);
		}
		return C;
	}

	friend Tensor operator-(const Tensor& A, const Tensor& B)
	{
		Tensor C(A.shape());
		for (int i = 0; i < A.shape().totalDimension(); i++)
		{
			C(i) = A(i) - B(i);
		}
		return C;
	}

	Tensor& operator+=(const Tensor& A);

	Tensor& operator-=(const Tensor& A);

	Tensor& operator*=(T a);

	Tensor& operator/=(T a);

	friend Tensor operator*(T a, const Tensor<T>& A) {
        Tensor B(A);
        B *= a;
        return B;
    }

	//////////////////////////////////////////////////////////
	// Adjust Dimensions
	//////////////////////////////////////////////////////////
	// Adjust Dimensions to a new TensorDim
	Tensor<T> AdjustDimensions(const TensorShape& newTDim)const;

	// Adjust the number of the active_ mode
	Tensor<T> AdjustActiveDim(size_t active, size_t mode)const;

	// Adjust the number of Tensors
	Tensor<T> AdjustStateDim(size_t n)const;

	// Reshape the tensor but keep the total size
	void Reshape(const TensorShape& tdim);

	//////////////////////////////////////////////////////////
	// Operations on Tensors
	//////////////////////////////////////////////////////////
	T singleDotProduct(const Tensor& A, size_t n, size_t m)const;

	//! A Dot Product between two Tensors.
	/*!
		This function calculates the doc-product between the tensors
		(*this) and A. The result is an overlap matrix sized with
		the number of tensors.
	 */

	Matrix<T> DotProduct(const Tensor<T>& A)const;

	/// This function will fill the Tensor with Zero-entries
	void Zero();

	// Getter for shape
	const TensorShape& shape()const { return shape_; }
	TensorShape& shape() { return shape_; }

protected:
	double conjugate(const double d) const {
		return d;
	}

	complex<double> conjugate(const complex<double> c) const {
		return conj(c);
	}

	TensorShape shape_;
	T* coeffs_;
	bool ownership_;
};

typedef Tensor<complex<double>> Tensorcd;
typedef Tensor<double> Tensord;

//////////////////////////////////////////////////////////
// Non-member functions
//////////////////////////////////////////////////////////
template<typename T>
T SingleDotProd(const Tensor<T>& A, const Tensor<T>& B, size_t n, size_t m);

template<typename T>
Tensor<T> productElementwise(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
void TensorContraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind);

template<typename T>
void Contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k, bool zero = true);

template<typename T>
Matrix<T> Contraction(const Tensor<T>& A, const Tensor<T>& B, size_t k);

template <typename T, typename U>
void MatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero = true);

template <typename T, typename U>
void TMatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero = true);

template <typename T, typename U>
void MatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template <typename T, typename U>
Tensor<T> MatrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

template <typename T, typename U>
void TensorMatrix(Tensor<T>& C, const Tensor<T>& B, const Matrix<U>& A, size_t mode, bool zero = true);

template<typename T, typename U>
Tensor<T> TensorMatrix(const Tensor<T>& B, const Matrix<U>& A, size_t mode);

template <typename T, typename U>
Tensor<T> multATB(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

template <typename T, typename U>
void multStateAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, bool zero = true);

template <typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B);

template<typename T, typename U>
void multStateArTB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B);

template <typename T, typename U>
Tensor<T> multStateArTB(const Matrix<U>& A, const Tensor<T>& B);

template <typename T, typename U>
void multAdd(Tensor<T>& A, const Tensor<T>& B, U coeff);

template<typename T>
void GramSchmidt(Tensor<T>& A);

//Projects B on A
template<typename T>
Tensor< complex<double> > Project(const Tensor< complex<double> >& A,
	const Tensor< T >& B);

template<typename T>
Tensor<T> ProjectOut(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Tensor< complex<double> > ProjectOrthogonal(const Tensor< complex<double> >& A,
	const Tensor< T >& B);

template<typename T>
Tensor<T> conj(Tensor<T> A);

template<typename T>
double Residual(Tensor<T> A, const Tensor<T>& B);

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A);

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B);

template<typename T>
ostream& operator<<(ostream& os, const Tensor<T>& A);

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A);

template<typename T>
bool operator==(const Tensor<T>& A, const Tensor<T>& B);

