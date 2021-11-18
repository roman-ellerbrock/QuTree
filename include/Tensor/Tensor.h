#pragma once
#include "TensorShape.h"

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
 * \brief This class represents an n-th order tensor
 *
 * The class provides tensor reshaping & flattening operations,
 * as well as specific cases of tensor-dot products that are
 * relevant in the context of tree tensor network states.
 * The tensor can be reinterpreted as a matrix by choosing a specific
 * target index k and by flattening all other indices into a compound
 * index I, resulting in A(k, I). The class uses this mapping to extend
 * matrix factorizations to tensors. If not other stated, the last index
 * is chosen and all previous indices are combined.
 *
 * Usage Examples:
 * TensorDim dim_({2, 3, 4, 1});
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

	explicit Tensor(const TensorShape& dim, Tensor<T>& ptr, bool ownership = true, bool InitZero = true);

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

	void write(ostream& os)const;

	void write(const string& filename)const;

	void read(istream& is);

	void read(const string& filename);

	//////////////////////////////////////////////////////////
	// Bracket Operators
	//////////////////////////////////////////////////////////
	T& operator()(size_t i)const;

	T& operator()(size_t i);

	const T& operator()(size_t i, size_t n)const;

	T& operator()(size_t i, size_t n);

	T& operator()(size_t bef, size_t j, size_t aft, size_t leaf);

	const T& operator()(size_t bef, size_t j, size_t aft, size_t leaf)const;

	T& operator()(const vector<size_t>& dims);

	const T& operator()(const vector<size_t>& dims) const;

	inline const T& operator[](const size_t idx)const {
		// Fast bracket operator
		return coeffs_[idx];
	}

	inline T& operator[](const size_t idx) {
		// Fast bracket operator
		return coeffs_[idx];
	}

	inline T& operator[](const vector<size_t>& idxs) {
		size_t I = indexMapping(idxs, shape_);
		return operator[](I);
	}

	//////////////////////////////////////////////////////////
	// Adjust Dimensions
	//////////////////////////////////////////////////////////
	// Adjust Dimensions to a new TensorDim
	Tensor<T> adjustDimensions(const TensorShape& newTDim)const;

	// Adjust the number of the active_ mode
	Tensor<T> adjustActiveDim(size_t active, size_t mode)const;

	// Adjust the number of Tensors
	Tensor<T> adjustStateDim(size_t n)const;

	// Reshape the tensor but keep the total size
	void reshape(const TensorShape& newShape);

	// Reshape the tensor and resize memory if required
	void resize(const TensorShape& newShape);

	//////////////////////////////////////////////////////////
	// Operations on Tensors
	//////////////////////////////////////////////////////////

	/// This function will fill the Tensor with Zero-entries
	void zero();

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
void elementwise(Tensor<T>& res, const Tensor<T>& A, const function<T(T)>& f);

template<typename T>
Tensor<T> elementwise(const Tensor<T>& A, const function<T(T)>& f);

template<typename T>
ostream& operator<<(ostream& os, const Tensor<T>& A);

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A);

template<typename T>
bool operator==(const Tensor<T>& A, const Tensor<T>& B);
