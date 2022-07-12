#pragma once
#include "TensorShape.h"
#include "rng.h"
#include <functional>
#include "hostMemory.h"

/**
 * \defgroup Core
 * \brief This group includes the basic datastructures in QuTree.
 *
 * These datastructures include the Vector, Matrix, and Tensor classes.
 */

template<typename T>
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
	Tensor()
		: coeffs_(new T[1]), ownership_(true) {}

	Tensor(const initializer_list<size_t>& dim, bool InitZero = true);

	// Constructor with TensorDim
	explicit Tensor(const TensorShape& dim, bool InitZero = true);

	// Construct from external memory
	explicit Tensor(const TensorShape& dim, T *ptr, bool ownership = true, bool InitZero = true);

	explicit Tensor(const TensorShape& dim, Tensor<T>& ptr, bool ownership = true, bool InitZero = true);

	explicit Tensor(istream& is);

	explicit Tensor(const string& filename);

	// Copy constructor
	Tensor(const Tensor& old);

	// Copy-Multiply constructor
	Tensor(const Tensor& old, T factor);

	// Move constructor
	Tensor(Tensor&& old) noexcept;

	// Copy Assignment Operator
	Tensor& operator=(const Tensor& old);

	// Move Assignment Operator
	Tensor& operator=(Tensor&& old) noexcept;

	// Destructor
	~Tensor();

	//////////////////////////////////////////////////////////
	// File handling
	//////////////////////////////////////////////////////////
	void print(ostream& os = cout) const;

	void write(ostream& os) const;

	void write(const string& filename) const;

	void read(istream& is);

	void read(const string& filename);

	//////////////////////////////////////////////////////////
	// Bracket Operators
	//////////////////////////////////////////////////////////
	T& operator()(size_t i) const;

	T& operator()(size_t i);

	const T& operator()(size_t i, size_t n) const;

	T& operator()(size_t i, size_t n);

	T& operator()(size_t bef, size_t j, size_t aft, size_t leaf);

	const T& operator()(size_t bef, size_t j, size_t aft, size_t leaf) const;

	T& operator()(const vector<size_t>& dims);

	const T& operator()(const vector<size_t>& dims) const;

	inline const T& operator[](const size_t idx) const {
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
	[[nodiscard]] Tensor<T> adjustDimensions(const TensorShape& newTDim) const;

	// Adjust the number of the active_ mode
	[[nodiscard]] Tensor<T> adjustActiveDim(size_t active, size_t mode) const;

	// Adjust the number of Tensors
	[[nodiscard]] Tensor<T> adjustStateDim(size_t n) const;

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
	T *data() { return coeffs_; }
	const T *data()const { return coeffs_; }
	polymorphic::hostMemory<T> host_{};
	bool ownership_;
private:
	T *coeffs_;
};

typedef Tensor<complex<double>> Tensorcd;

typedef Tensor<double> Tensord;

//////////////////////////////////////////////////////////
// Non-member functions
//////////////////////////////////////////////////////////

template<typename T>
void elementwise(Tensor<T>& res, const Tensor<T>& A, const std::function<T(T)>& f);

template<typename T>
Tensor<T> elementwise(const Tensor<T>& A, const std::function<T(T)>& f);

template<typename T>
ostream& operator<<(ostream& os, const Tensor<T>& A);

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A);

template<typename T>
bool operator==(const Tensor<T>& A, const Tensor<T>& B);


template<typename T>
[[nodiscard]] Tensor<T> random(const TensorShape& shape, mt19937& gen);

auto randomcd = [](const TensorShape& shape, mt19937& gen = rng::gen) {
	return random<complex<double>>(shape, gen);
};
auto randomd = [](const TensorShape& shape, mt19937& gen = rng::gen) {
	return random<double>(shape, gen);
};

template<typename T>
[[nodiscard]] Tensor<T> randomGen(const TensorShape& shape);

auto randomGencd = [](const TensorShape& shape) {
	return random<complex<double>>(shape, rng::gen);
};
auto randomGend = [](const TensorShape& shape) {
	return random<double>(shape, rng::gen);
};

template<typename T>
[[nodiscard]] Tensor<T> arange(const TensorShape& shape);

constexpr auto arangecd = arange<complex<double>>;
constexpr auto aranged = arange<double>;


template<typename T>
[[nodiscard]] Tensor<T> identity(const TensorShape& shape);

constexpr auto identitycd = identity<complex<double>>;
constexpr auto identityd = identity<double>;


template<typename T>
[[nodiscard]] Tensor<T> delta(const TensorShape& shape);

constexpr auto deltacd = delta<complex<double>>;
constexpr auto deltad = delta<double>;


template <typename T>
using Matrix = Tensor<T>;

typedef Matrix<complex<double>> Matrixcd;

typedef Matrix<double> Matrixd;


template <typename T>
using Vector = Tensor<T>;

typedef Vector<complex<double>> Vectorcd;

typedef Vector<double> Vectord;

