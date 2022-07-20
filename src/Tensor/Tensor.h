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

template<typename T, template <typename> class Memory = polymorphic::hostMemory>
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
	Tensor() = default;

	Tensor(const initializer_list<size_t>& dim, bool InitZero = true);

	// Constructor with TensorDim
	explicit Tensor(const TensorShape& dim, bool InitZero = true);

	explicit Tensor(istream& is);

	explicit Tensor(const string& filename);
	
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
	const T& operator()(size_t i) const;

	T& operator()(size_t i);

	const T& operator()(size_t i, size_t n) const;

	T& operator()(size_t i, size_t n);

	T& operator()(size_t bef, size_t j, size_t aft, size_t leaf);

	const T& operator()(size_t bef, size_t j, size_t aft, size_t leaf) const;

	T& operator()(const vector<size_t>& dims);

	const T& operator()(const vector<size_t>& dims) const;

	inline const T& operator[](const size_t idx) const {
		// Fast bracket operator
		return data()[idx];
	}

	inline T& operator[](const size_t idx) {
		// Fast bracket operator
		return data()[idx];
	}

	inline T& operator[](const vector<size_t>& idxs) {
		size_t I = indexMapping(idxs, shape_);
		return operator[](I);
	}

	//////////////////////////////////////////////////////////
	// Adjust Dimensions
	//////////////////////////////////////////////////////////
	// Adjust Dimensions to a new TensorDim
	[[nodiscard]] Tensor<T, Memory> adjustDimensions(const TensorShape& newTDim) const;

	// Adjust the number of the active_ mode
	[[nodiscard]] Tensor<T, Memory> adjustActiveDim(size_t active, size_t mode) const;

	// Adjust the number of Tensors
	[[nodiscard]] Tensor<T, Memory> adjustStateDim(size_t n) const;

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
	auto *data() { return mem_.data(); }
	const auto*data()const { return mem_.data(); }
	Memory<T>& mem() { return mem_; }
	const Memory<T>& mem() const { return mem_; }
private:
	Memory<T> mem_{};
//	T *coeffs_;
};

typedef Tensor<complex<double>> Tensorcd;

typedef Tensor<double> Tensord;

//////////////////////////////////////////////////////////
// Non-member functions
//////////////////////////////////////////////////////////

template<typename T, template <typename> class Mem>
void elementwise(Tensor<T,Mem>& res, const Tensor<T,Mem>& A, const std::function<T(T)>& f);

template<typename T, template <typename> class Mem>
Tensor<T,Mem> elementwise(const Tensor<T,Mem>& A, const std::function<T(T)>& f);

template<typename T, template <typename> class Mem>
ostream& operator<<(ostream& os, const Tensor<T,Mem>& A);

template<typename T, template <typename> class Mem>
istream& operator>>(istream& is, Tensor<T,Mem>& A);

template<typename T, template <typename> class Mem>
bool operator==(const Tensor<T,Mem>& A, const Tensor<T,Mem>& B);


template<typename T, template <typename> class Mem = polymorphic::hostMemory>
[[nodiscard]] Tensor<T,Mem> random(const TensorShape& shape, mt19937& gen);

auto randomcd = [](const TensorShape& shape, mt19937& gen = rng::gen) {
	return random<complex<double>>(shape, gen);
};
auto randomd = [](const TensorShape& shape, mt19937& gen = rng::gen) {
	return random<double>(shape, gen);
};

template<typename T, template <typename> class Mem>
[[nodiscard]] Tensor<T,Mem> randomGen(const TensorShape& shape);

auto randomGencd = [](const TensorShape& shape) {
	return random<complex<double>>(shape, rng::gen);
};
auto randomGend = [](const TensorShape& shape) {
	return random<double>(shape, rng::gen);
};

template<typename T, template <typename> class Mem = polymorphic::hostMemory>
[[nodiscard]] Tensor<T, Mem> arange(const TensorShape& shape);

constexpr auto arangecd = arange<complex<double>>;
constexpr auto aranged = arange<double>;


template<typename T, template <typename> class Mem = polymorphic::hostMemory>
[[nodiscard]] Tensor<T,Mem> identity(const TensorShape& shape);

constexpr auto identitycd = identity<complex<double>>;
constexpr auto identityd = identity<double>;


template<typename T, template <typename> class Mem = polymorphic::hostMemory>
[[nodiscard]] Tensor<T,Mem> delta(const TensorShape& shape);

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
