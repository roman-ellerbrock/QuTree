#pragma once
#include "TensorDim.h"
#include "Core/Matrix.h"


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
	Tensor() :coeffs(new T[1]) {}

	// Constructor with TensorDim
	explicit Tensor(const TensorDim& dim_, bool InitZero = true);

	explicit Tensor(istream& is);

	explicit Tensor(const string& filename);

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

	T& operator()(size_t i, size_t n)const;

	T& operator()(size_t i, size_t n);

	T& operator()(size_t i, size_t j, size_t k, size_t f, size_t n);

	T& operator()(size_t i, size_t j, size_t k, size_t f, size_t n)const;

	// double hole operator
	T& operator()(size_t bef, size_t i, size_t mid, size_t j, size_t beh,
		size_t mode1, size_t mode2, size_t n)const;

	inline T& operator[](const size_t idx)const
	{
		// Fast bracket operator
		return coeffs[idx];
	}

	inline T& operator[](const size_t idx)
	{
		// Fast bracket operator
		return coeffs[idx];
	}

	//////////////////////////////////////////////////////////
	// Math Operators
	//////////////////////////////////////////////////////////
	friend Tensor<T> operator+(const Tensor<T>& A, const Tensor<T>& B)
	{
		assert(A.Dim().GetDimTot() == B.Dim().GetDimTot());
		Tensor C(A.Dim());
		for (int i = 0; i < A.Dim().GetDimTot(); i++)
		{
			C(i) = A(i) + B(i);
		}
		return C;
	}

	friend Tensor operator-(const Tensor& A, const Tensor& B)
	{
		Tensor C(A.Dim());
		for (int i = 0; i < A.Dim().GetDimTot(); i++)
		{
			C(i) = A(i) - B(i);
		}
		return C;
	}

	void operator+=(const Tensor& A);

	void operator-=(const Tensor& A);

	void operator*=(T a);

    friend Tensor operator*(T a, const Tensor<T>& A) {
        Tensor B(A);
        B *= a;
        return B;
    }
    void operator/=(T a);

	Tensor coeffprod(const Tensor& A, const Tensor& B);

	//////////////////////////////////////////////////////////
	// Adjust Dimensions
	//////////////////////////////////////////////////////////
	// Adjust Dimensions to a new TensorDim
	Tensor<T> AdjustDimensions(const TensorDim& newTDim)const;

	// Adjust the number of the active_ mode
	Tensor<T> AdjustActiveDim(size_t active, size_t mode)const;

	// Adjust the number of Tensors
	Tensor<T> AdjustStateDim(size_t n)const;

	// Reshape the tensor but keep the total size
	void Reshape(const TensorDim& tdim);

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

	// Getter for Dim
	const TensorDim& Dim()const { return dim; }
	TensorDim& Dim() { return dim; }

protected:
	double conjugate(const double d) const {
		return d;
	}

	complex<double> conjugate(const complex<double> c) const {
		return conj(c);
	}

	TensorDim dim;
	T* coeffs;
};

typedef Tensor<complex<double>> Tensorcd;
typedef Tensor<double> Tensord;

//////////////////////////////////////////////////////////
// Non-member functions
//////////////////////////////////////////////////////////
template<typename T>
T SingleDotProd(const Tensor<T>& A, const Tensor<T>& B, size_t n, size_t m);

template<typename T>
void TensorHoleProduct(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind);

/*template<typename T>
void HoleProduct(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k);

template<typename T>
Matrix<T> HoleProduct(const Tensor<T>& A, const Tensor<T>& B, size_t k);
 */

template <typename T, typename U>
void MatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero = true);

template <typename T, typename U>
void Tmattensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind, bool zero = true);

template <typename T, typename U>
void multAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template <typename T, typename U>
Tensor<T> multAB(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

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
ostream& operator<<(ostream& os, const Tensor<T>& A);

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A);

template<typename T>
bool operator==(const Tensor<T>& A, const Tensor<T>& B);

