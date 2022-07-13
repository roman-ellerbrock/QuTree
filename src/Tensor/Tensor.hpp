#pragma once
#include "Tensor.h"
#include "stdafx.h"
#include <cstring>

template<typename T, template <typename> class Mem>
Tensor<T,Mem>::Tensor(const initializer_list<size_t>& dims, bool InitZero)
	:Tensor(TensorShape(dims), InitZero) {}

template<typename T, template <typename> class Mem>
Tensor<T,Mem>::Tensor(const TensorShape& shape, const bool InitZero)
	:shape_(shape), mem_(shape.totalDimension()) {
	if (InitZero) { zero(); }
}

template<typename T, template <typename> class Mem>
Tensor<T,Mem>::Tensor(istream& is)
	:Tensor() {
	read(is);
}

template<typename T, template <typename> class Mem>
Tensor<T,Mem>::Tensor(const string& filename)
	: Tensor() {
	ifstream is(filename);
	read(is);
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::zero() {
	memset(data(), 0, shape_.totalDimension() * sizeof(T));
}

//////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////

template<typename T, template <typename> class Mem>
inline const T& Tensor<T,Mem>::operator()(const size_t i) const {
	return data()[i];
}

template<typename T, template <typename> class Mem>
inline T& Tensor<T,Mem>::operator()(const size_t i) {
	return data()[i];
}

//////////////////////////////////////////////////////////
// Bracket Operators
//////////////////////////////////////////////////////////
template<typename T, template <typename> class Mem>
inline const T& Tensor<T,Mem>::operator()(const size_t i, const size_t n) const {
	size_t dimpart = shape_.lastBefore();
	return data()[n * dimpart + i];
}

template<typename T, template <typename> class Mem>
inline T& Tensor<T,Mem>::operator()(const size_t i, const size_t n) {
	size_t dimpart = shape_.lastBefore();
	return data()[n * dimpart + i];
}

template<typename T, template <typename> class Mem>
inline T& Tensor<T,Mem>::operator()(size_t bef, size_t i, size_t aft, size_t leaf) {
	size_t before = shape_.before(leaf);
	size_t dim = shape_[leaf];
	size_t idx = aft * before * dim + i * before + bef;
	return data()[idx];
}

template<typename T, template <typename> class Mem>
inline const T& Tensor<T,Mem>::operator()(size_t bef, size_t i, size_t aft, size_t leaf) const {
	size_t before = shape_.before(leaf);
	size_t dim = shape_[leaf];
	size_t idx = aft * before * dim + i * before + bef;
	return data()[idx];
}

template<typename T, template <typename> class Mem>
T& Tensor<T,Mem>::operator()(const vector<size_t>& dims) {
	return operator()(indexMapping(dims, shape_));
}

template<typename T, template <typename> class Mem>
const T& Tensor<T,Mem>::operator()(const vector<size_t>& dims) const {
	return operator()(indexMapping(dims, shape_));
}

//////////////////////////////////////////////////////////
// File handling
//////////////////////////////////////////////////////////
template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::print(ostream& os) const {
	if (shape_.empty()) {
		os << "[ ]" << endl;
		return;
	}
	for (size_t n = 0; n < shape_.lastDimension(); n++) {
		for (size_t i = 0; i < shape_.lastBefore(); i++)
			os << (*this)(i, n) << " ";
		os << endl;
	}
	os << endl;
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::write(ostream& os) const {
	// Verification
	os.write("TENS", 4);

	// Write the TensorDim
	shape_.write(os);

	// Write the size
	int32_t size = sizeof(T);
	os.write((char *) &size, sizeof(size));

	// Write the Coefficients
	for (size_t i = 0; i < shape_.totalDimension(); i++) {
		T Coeff_now = operator()(i);
		os.write((char *) &Coeff_now, size);
	}
	os.flush();
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::write(const string& file) const {
	ofstream os(file);
	write(os);
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::read(istream& is) {
	// Check if binary string contains a Tensor
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TENS");
	assert(s_key == s_check);

	// Read the TensorDim
	TensorShape newtdim;
	newtdim.readDim(is);

	// Resize the Tensor
	(*this) = Tensor<T,Mem>(newtdim, false);

	// Read the size
	int32_t size;
	is.read((char *) &size, sizeof(size));
	assert(size == sizeof(T));

	// Read the coefficients

	for (size_t i = 0; i < shape_.totalDimension(); i++) {
		T Coeff_now;
		is.read((char *) &Coeff_now, size);
		operator()(i) = Coeff_now;
	}
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

//////////////////////////////////////////////////////////
// Adjust Dimensions
//////////////////////////////////////////////////////////
template<typename T, template <typename> class Mem>
Tensor<T,Mem> Tensor<T,Mem>::adjustDimensions(const TensorShape& newTDim) const {
	// Increase the dimensions of the Tensor from old TensorDim
	// to new TensorDim 

	assert(newTDim.order() == shape_.order());
	// Increase the active_ modes
	Tensor<T,Mem> Acoeff(*this);
	for (size_t k = 0; k < shape_.order(); k++) {
		size_t act = newTDim[k];
		Acoeff = Acoeff.adjustActiveDim(act, k);
	}

	// Increase the number of Tensors
	size_t ntens = newTDim.lastDimension();
	Acoeff = Acoeff.adjustStateDim(ntens);

	return Acoeff;
}

template<typename T, template <typename> class Mem>
Tensor<T,Mem> Tensor<T,Mem>::adjustActiveDim(size_t active, size_t mode) const {
	// Adjust the active_ dimension in the coordinate "mode".
	// If the new active_ is smaller, the norm of the tensors is
	// not conserved.

	assert(mode < shape_.order());

	// Create a new Tensor with the adjusted dim_
	vector<size_t> dimlist = shape_.dimensions();
	dimlist[mode] = active;
	TensorShape newTDim(dimlist);
	Tensor<T,Mem> newT(newTDim);

	// Copy the coefficients
	size_t before = shape_.before(mode);
	size_t after = shape_.after(mode);
	size_t minactive = min(active, shape_[mode]);
	/// Offsets are used to add new & delete functions at first indices.
	/// This ensures low-to-high occupancy convention.
	size_t offset_old = shape_[mode] - minactive;
	size_t offset_new = active - minactive;
	for (size_t l = 0; l < after; l++) {
		for (size_t j = 0; j < minactive; j++) {
			for (size_t i = 0; i < before; i++) {
				newT(i, j + offset_new, l, mode) = operator()(i, j + offset_old, l, mode);
			}
		}
	}

	return newT;
}

// Adjust the size of Tensor 
template<typename T, template <typename> class Mem>
Tensor<T,Mem> Tensor<T,Mem>::adjustStateDim(size_t n) const {
	return adjustActiveDim(n, shape_.lastIdx());
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::reshape(const TensorShape& newShape) {
	/// Check that total size is the same
	assert(shape_.totalDimension() == newShape.totalDimension());
	shape_ = newShape;
}

template<typename T, template <typename> class Mem>
void Tensor<T,Mem>::resize(const TensorShape& newShape) {
	/// resize if required
	if (shape_.totalDimension() != newShape.totalDimension()) {
		mem_ = polymorphic::hostMemory<T>(newShape.totalDimension());
	}
	shape_ = newShape;
}

//////////////////////////////////////////////////////////
// Operations on Tensors
//////////////////////////////////////////////////////////
///  f(A(i))
template<typename T, template <typename> class Mem>
void elementwise(Tensor<T,Mem>& res, const Tensor<T,Mem>& A, const function<T(T)>& f) {
	assert(A.Dim1() == res.Dim1());
	assert(A.Dim2() == res.Dim2());
	for (size_t i = 0; i < A.Dim1() * A.Dim2(); ++i) {
		res[i] = f(A[i]);
	}
}

///  = f(A(i))
template<typename T, template <typename> class Mem>
Tensor<T,Mem> elementwise(const Tensor<T,Mem>& A, const function<T(T)>& f) {
	Tensor<T,Mem> res(A.Dim1(), A.Dim2(), false);
	elementwise(res, A, f);
	return res;
}

//////////////////////////////////////////////////////////
/// Non-member functions
//////////////////////////////////////////////////////////
template<typename T, template <typename> class Mem>
ostream& operator<<(ostream& os, const Tensor<T,Mem>& A) {
	A.write(os);
	return os;
}

template<typename T, template <typename> class Mem>
istream& operator>>(istream& is, Tensor<T,Mem>& A) {
	A.read(is);
	return is;
}

template<typename T, template <typename> class Mem>
bool operator==(const Tensor<T,Mem>& A, const Tensor<T,Mem>& B) {
	if (A.shape_ != B.shape_) { return false; }
	for (size_t k = 0; k < A.shape_.totalDimension(); ++k) {
		if (A[k] != B[k]) { return false; }
	}
	return true;
}

template<typename T, template <typename> class Mem>
Tensor<T,Mem> random(const TensorShape& shape, mt19937& gen) {
	Tensor<T,Mem> A(shape, false);
	normal_distribution dist(0., 1.);
	for (size_t i = 0; i < shape.totalDimension(); ++i) {
		A[i] = dist(gen);
	}
	return A;
}

template<typename T, template <typename> class Mem>
Tensor<T,Mem> randomGen(const TensorShape& shape) {
	return random<T>(shape, rng::gen);
}

template<typename T, template <typename> class Mem>
[[nodiscard]] Tensor<T,Mem> arange(const TensorShape& shape) {
	Tensor<T,Mem> A(shape, false);
	for (size_t i = 0; i < shape.totalDimension(); ++i) {
		A[i] = i;
	}
	return A;
}

template<typename T, template <typename> class Mem>
[[nodiscard]] Tensor<T,Mem> identity(const TensorShape& shape) {
	Tensor<T,Mem> A(shape, true);
	size_t n = (shape.lastDimension() < shape.lastBefore()) ? shape.lastDimension() : shape.lastBefore();
	for (size_t i = 0; i < n; ++i) {
		A(i, i) = 1.;
	}
	return A;
}

template<typename T, template <typename> class Mem>
[[nodiscard]] Tensor<T,Mem> delta(const TensorShape& shape) {
	Tensor<T,Mem> A(shape);
	A(0) = 1.;
	return A;
}

