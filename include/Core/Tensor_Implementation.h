#pragma once
#include "Tensor.h"
#include "TensorShape.h"
#include "stdafx.h"
//#include <omp.h> //TODO: have this here by default?

template<typename T>
Tensor<T>::Tensor(const initializer_list<size_t>& dims, bool InitZero)
	:Tensor(TensorShape(dims), InitZero) {}

template<typename T>
Tensor<T>::Tensor(const TensorShape& dim, T *ptr, bool ownership, bool InitZero)
	:shape_(dim), coeffs_(ptr), ownership_(ownership) {
	if (InitZero) { zero(); }
}

// Construct from external tensor that holds the memory
template<typename T>
Tensor<T>::Tensor(const TensorShape& dim, Tensor<T>& A, bool ownership, bool InitZero)
	: Tensor<T>(dim, &A[0], ownership, false){
	if (dim.totalDimension() > A.shape_.totalDimension()) {
		cerr << "Error: memory too small for tensor.\n";
		ownership_ = false;
		exit(1);
	}
	if (ownership &! A.ownership_) {
		cerr << "Error: cannot transfer ownership, since original tensor was not owning memory.\n";
		ownership_ = false;
		exit(1);
	}
	if (ownership) {
		A.ownership_ = false;
	}
	if (InitZero) { zero(); }
}


template<typename T>
Tensor<T>::Tensor(const TensorShape& dim, const bool InitZero)
	:shape_(dim), coeffs_((T *) malloc(dim.totalDimension() * sizeof(T))), ownership_(true) {
//	:shape_(dim), coeffs_(new T[dim.totalDimension()]), ownership_(true) {
	if (InitZero) { zero(); }
}

template<typename T>
Tensor<T>::Tensor(istream& is)
	:Tensor() {
	read(is);
}

template<typename T>
Tensor<T>::Tensor(const string& filename)
	: Tensor() {
	ifstream is(filename);
	read(is);
}

// Copy constructor
template<typename T>
Tensor<T>::Tensor(const Tensor& old)
	:Tensor(old.shape_, false) {
	memcpy(coeffs_, old.coeffs_, shape_.totalDimension() * sizeof(T));
}

// Copy-Multyply constructor
template<typename T>
Tensor<T>::Tensor(const Tensor& old, T factor)
	:Tensor(old.shape_, false) {
	for (size_t i = 0; i < shape_.totalDimension(); i++) {
		coeffs_[i] = old.coeffs_[i] * factor;
	}
}

// Move constructor
template<typename T>
Tensor<T>::Tensor(Tensor&& old) noexcept
	:shape_(old.shape_), coeffs_(old.coeffs_), ownership_(old.ownership_) {
	old.coeffs_ = nullptr;
	old.ownership_ = false;
}

// Copy Assignment Operator
template<typename T>
Tensor<T>& Tensor<T>::operator=(const Tensor& old) {
	if (this == &old) {
		return *this;
	} else if (old.shape_ == this->shape_) {
		memcpy(coeffs_, old.coeffs_, shape_.totalDimension() * sizeof(T));
	} else {
		Tensor tmp(old);
		*this = move(tmp);
	}
	return *this;
}

// Move Assignment Operator
template<typename T>
Tensor<T>& Tensor<T>::operator=(Tensor&& old) noexcept {
	delete[] coeffs_;
	shape_ = old.shape_;
	coeffs_ = old.coeffs_;
	ownership_ = old.ownership_;
	old.coeffs_ = nullptr;
	old.ownership_ = false;
	return *this;
}

template<typename T>
Tensor<T>::~Tensor() {
	if (ownership_) { delete[] coeffs_; }
}

//////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////

template<typename T>
inline T& Tensor<T>::operator()(const size_t i) const {
	size_t dimtot = shape_.totalDimension();
	assert(i < dimtot);
	return coeffs_[i];
}

template<typename T>
inline T& Tensor<T>::operator()(const size_t i) {
	size_t dimtot = shape_.totalDimension();
	assert(i < dimtot);
	return coeffs_[i];
}

//////////////////////////////////////////////////////////
// Bracket Operators
//////////////////////////////////////////////////////////
template<typename T>
inline const T& Tensor<T>::operator()(const size_t i, const size_t n) const {
	size_t dimpart = shape_.lastBefore();
	assert(i < dimpart);
	assert(n < shape_.lastDimension());
	return coeffs_[n * dimpart + i];
}

template<typename T>
inline T& Tensor<T>::operator()(const size_t i, const size_t n) {
	size_t dimpart = shape_.lastBefore();
	assert(i < dimpart);
	assert(n < shape_.lastDimension());
	return coeffs_[n * dimpart + i];
}

template<typename T>
inline T& Tensor<T>::operator()(size_t bef, size_t i, size_t aft, size_t leaf) {
	assert(leaf < shape_.order());
	assert(bef < shape_.before(leaf));
	assert(i < shape_[leaf]);
	assert(aft < shape_.after(leaf));
	size_t before = shape_.before(leaf);
	size_t dim = shape_[leaf];
	size_t idx = aft * before * dim + i * before + bef;
	// @TODO: remove when tested
	assert(idx < shape_.totalDimension());
	return coeffs_[idx];
}

template<typename T>
inline const T& Tensor<T>::operator()(size_t bef, size_t i, size_t aft, size_t leaf) const {
	assert(leaf < shape_.order());
	assert(bef < shape_.before(leaf));
	assert(i < shape_[leaf]);
	assert(aft < shape_.after(leaf));
	size_t before = shape_.before(leaf);
	size_t dim = shape_[leaf];
	size_t idx = aft * before * dim + i * before + bef;
	// @TODO: remove when tested
	assert(idx < shape_.totalDimension());
	return coeffs_[idx];
}

template<typename T>
T& Tensor<T>::operator()(const vector<size_t>& dims) {
	return operator()(indexMapping(dims, shape_));
}

template<typename T>
const T& Tensor<T>::operator()(const vector<size_t>& dims) const {
	return operator()(indexMapping(dims, shape_));
}

//////////////////////////////////////////////////////////
// File handling
//////////////////////////////////////////////////////////
template<typename T>
void Tensor<T>::print(ostream& os) const {
	for (size_t n = 0; n < shape_.lastDimension(); n++) {
		for (size_t i = 0; i < shape_.lastBefore(); i++)
			os << (*this)(i, n) << " ";
		os << endl;
	}
	os << endl;
}

template<typename T>
void Tensor<T>::write(ostream& os) const {
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

template<typename T>
void Tensor<T>::write(const string& file) const {
	ofstream os(file);
	write(os);
}

template<typename T>
void Tensor<T>::read(istream& is) {
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
	(*this) = Tensor<T>(newtdim, false);

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

template<typename T>
void Tensor<T>::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

//////////////////////////////////////////////////////////
// Math Operators
//////////////////////////////////////////////////////////

template<typename T>
Tensor<T>& Tensor<T>::operator+=(const Tensor& A) {
	assert(A.shape_.totalDimension() == shape_.totalDimension());
	T const *Ax = A.coeffs_;
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		coeffs_[i] += Ax[i];
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator-=(const Tensor& A) {
	assert(A.shape_.totalDimension() == shape_.totalDimension());
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		(*this)(i) -= A(i);
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator*=(T a) {
	for (size_t i = 0; i < shape_.totalDimension(); i++) {
		operator()(i) = a * operator()(i);
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator/=(T a) {
	for (size_t i = 0; i < shape_.totalDimension(); i++) {
		operator()(i) = operator()(i) / a;
	}
	return *this;
}

template<typename T>
Tensor<T> productElementwise(const Tensor<T>& A, const Tensor<T>& B) {
	assert(A.shape_.totalDimension() == B.shape_.totalDimension());
	Tensor<T> C(A.shape_);
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		C(i) = A(i) * B(i);
	}
	return C;
}

//////////////////////////////////////////////////////////
// Adjust Dimensions
//////////////////////////////////////////////////////////
template<typename T>
Tensor<T> Tensor<T>::adjustDimensions(const TensorShape& newTDim) const {
	// Increase the dimensions of the Tensor from old TensorDim
	// to new TensorDim 

	assert(newTDim.order() == shape_.order());
	// Increase the active_ modes
	Tensor<T> Acoeff(*this);
	for (size_t k = 0; k < shape_.order(); k++) {
		size_t act = newTDim[k];
		Acoeff = Acoeff.adjustActiveDim(act, k);
	}

	// Increase the number of Tensors
	size_t ntens = newTDim.lastDimension();
	Acoeff = Acoeff.adjustStateDim(ntens);

	return Acoeff;
}

template<typename T>
Tensor<T> Tensor<T>::adjustActiveDim(size_t active, size_t mode) const {
	// Adjust the active_ dimension in the coordinate "mode".
	// If the new active_ is smaller, the norm of the tensors is
	// not conserved.

	assert(mode < shape_.order());

	// Create a new Tensor with the adjusted dim_
	vector<size_t> dimlist = shape_.dimensions();
	dimlist[mode] = active;
	TensorShape newTDim(dimlist);
	Tensor<T> newT(newTDim);

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
template<typename T>
Tensor<T> Tensor<T>::adjustStateDim(size_t n) const {
	return adjustActiveDim(n, shape_.lastIdx());
}

template<typename T>
void Tensor<T>::reshape(const TensorShape& new_dim) {
	/// Check that total size is the same
	assert(shape_.totalDimension() == new_dim.totalDimension());
	shape_ = new_dim;
}

//////////////////////////////////////////////////////////
// Operations on Tensors
//////////////////////////////////////////////////////////
template<typename T>
void Tensor<T>::zero() {
	memset(coeffs_, 0, shape_.totalDimension() * sizeof(T));
}

//////////////////////////////////////////////////////////
/// Non-member functions
//////////////////////////////////////////////////////////
template<typename T>
T singleDotProd(const Tensor<T>& A, const Tensor<T>& B, size_t n, size_t m) {
	TensorShape tdima(A.shape_);
	TensorShape tdimb(B.shape_);

	size_t nmax = tdima.lastDimension();
	size_t mmax = tdimb.lastDimension();
	size_t npart = tdima.lastBefore();

	// Every tensor can have different amount of states but same dimpart
	assert(npart == tdimb.lastBefore());
	assert(n < nmax);
	assert(m < mmax);

	T result = 0;
#pragma omp parallel for reduction(+:result)
	for (size_t i = 0; i < npart; i++) {
		result += conj(A(i, n)) * B(i, m);
	}
	return result;
}


template<typename T, typename U>
void multAdd(Tensor<T>& A, const Tensor<T>& B, U coeff) {
	const TensorShape& tdim = A.shape_;
	const TensorShape& tdim_2 = A.shape_;
	size_t dimtot = tdim.totalDimension();
	assert(dimtot == tdim_2.totalDimension());
	for (size_t i = 0; i < dimtot; ++i) {
		A(i) += coeff * B(i);
	}
}

template<typename T>
void gramSchmidt(Tensor<T>& A) {
	// @TODO: Fill in auto-refill

	// control parameters
	size_t maxiter = 15;
	double conver = 1e-12;
	double errorconver = 1e-9;

	TensorShape tdim(A.shape_);
	size_t ntensor = tdim.lastDimension();
	size_t dimpart = tdim.lastBefore();

	for (size_t n = 0; n < ntensor; n++) {
		size_t iter = 0;
		double accumoverlap = 1.;
		// orthogonalize on all previous ones and then normalize
		while ((accumoverlap > conver) && (iter < maxiter)) {
			iter++;
			accumoverlap = 0;
			for (size_t m = 0; m < n; m++) {
				// orthogonalize
				T overlap = singleDotProd(A, A, m, n);
				accumoverlap += abs(overlap);
				for (size_t i = 0; i < dimpart; i++) {
					A(i, n) -= overlap * A(i, m);
				}
			}

			// renormalize
			T norm = singleDotProd(A, A, n, n);
			if (abs(norm) != 0) {
				norm = sqrt(real(norm));
				for (size_t i = 0; i < dimpart; i++) {
					A(i, n) /= norm;
				}
			}
		}
		// Error message
		if (accumoverlap >= errorconver) {
			cout << "Error: No orthogonality in Gram-Schmidt" << endl;
			cout << "Error measurement: " << conver << endl;
			cout << "Present error: " << accumoverlap << endl;
			cout << "Error acceptance: " << errorconver << endl;

			assert(0);
		}
	}
}

template<typename T>
void gramSchmidt(Tensor<T>& A, size_t k) {

}

template<typename T>
Tensor<T> conj(Tensor<T> A) {
	return elementwise(A, conj);
}

template<typename T>
double residual(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	return A.norm();
}

template<typename T>
ostream& operator<<(ostream& os, const Tensor<T>& A) {
	A.write(os);
	return os;
}

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A) {
	A.read(is);
	return is;
}

template<typename T>
bool operator==(const Tensor<T>& A, const Tensor<T>& B) {
	if (A.shape_ != B.shape_) { return false; }
	for (size_t k = 0; k < A.shape_.totalDimension(); ++k) {
		if (A[k] != B[k]) { return false; }
	}
	return true;
}

template<typename T>
void elementwise(Tensor<T>& res, const Tensor<T>& A, const function<T(T)>& f) {
	assert(A.Dim1() == res.Dim1());
	assert(A.Dim2() == res.Dim2());
	for (size_t i = 0; i < A.Dim1() * A.Dim2(); ++i) {
		res[i] = f(A[i]);
	}
}

template<typename T>
Tensor<T> elementwise(const Tensor<T>& A, const function<T(T)>& f) {
	Tensor<T> res(A.Dim1(), A.Dim2(), false);
	elementwise(res, A, f);
	return res;
}
