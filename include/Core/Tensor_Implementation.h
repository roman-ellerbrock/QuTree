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
	if (dim.totalDimension() > A.shape().totalDimension()) {
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
	memcpy(coeffs_, old.coeffs_, shape().totalDimension() * sizeof(T));
//	for (size_t i = 0; i < shape_.totalDimension(); i++) {
//		coeffs_[i] = old.coeffs_[i];
//	}
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
	} else if (old.shape() == this->shape()) {
		memcpy(coeffs_, old.coeffs_, shape().totalDimension() * sizeof(T));
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
	assert(A.shape().totalDimension() == shape().totalDimension());
	T const *Ax = A.coeffs_;
	for (size_t i = 0; i < A.shape().totalDimension(); i++) {
		coeffs_[i] += Ax[i];
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator-=(const Tensor& A) {
	assert(A.shape().totalDimension() == shape().totalDimension());
	for (size_t i = 0; i < A.shape().totalDimension(); i++) {
		(*this)(i) -= A(i);
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator*=(T a) {
	for (size_t i = 0; i < shape().totalDimension(); i++) {
		operator()(i) = a * operator()(i);
	}
	return *this;
}

template<typename T>
Tensor<T>& Tensor<T>::operator/=(T a) {
	for (size_t i = 0; i < shape().totalDimension(); i++) {
		operator()(i) = operator()(i) / a;
	}
	return *this;
}

template<typename T>
Tensor<T> productElementwise(const Tensor<T>& A, const Tensor<T>& B) {
	assert(A.shape().totalDimension() == B.shape().totalDimension());
	Tensor<T> C(A.shape());
	for (size_t i = 0; i < A.shape().totalDimension(); i++) {
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
	return adjustActiveDim(n, shape().lastIdx());
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
Matrix<T> Tensor<T>::dotProduct(const Tensor<T>& A) const {
	TensorShape tdima(A.shape());
	size_t nmax = tdima.lastDimension();
	size_t mmax = shape_.lastDimension();
	Matrix<T> S(mmax, nmax);
	contraction(S, *this, A, shape_.lastIdx());
	return S;
}

template<typename T>
void Tensor<T>::zero() {
	memset(coeffs_, 0, shape_.totalDimension() * sizeof(T));
}

template<typename T>
Tensor<T>::Tensor(const Matrix<T>& mat)
	: Tensor<T>({mat.dim1(), mat.dim2()}) {

	for (size_t i = 0; i < mat.dim2(); ++i) {
		for (size_t k = 0; k < mat.dim1(); ++k) {
			this->operator[](indexMapping({k, i}, shape_)) = mat(k, i);
		}
	}
}

//////////////////////////////////////////////////////////
/// Non-member functions
//////////////////////////////////////////////////////////
template<typename T>
T singleDotProd(const Tensor<T>& A, const Tensor<T>& B, size_t n, size_t m) {
	TensorShape tdima(A.shape());
	TensorShape tdimb(B.shape());

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

extern "C" {
// subroutine matvec (mulpsi, psi, matrix, a, b, c, add)
// subroutine rhomat (bra,ket,matrix,a,b,c)
void matvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void ctmatvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void rmatvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void rhomat_(double *Bra, double *Ket, double *M,
	int *a, int *b, int *c);
}

template<typename T, typename U>
void contraction1(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {

	if (zero) { h.zero(); }

	for (size_t a = 0; a < A; ++a) {
		for (size_t b = 0; b < B; ++b) {
			for (size_t b2 = 0; b2 < B2; ++b2) {
				for (size_t c = 0; c < C; ++c) {
//					h(b, b2) += conj(bra(a, b, c)) * ket(a, b2, c);
					h(b, b2) += conj(bra[a + b * A + c * A * B]) * ket(a + b2 * A + c * A * B2);
				}
			}
		}
	}
}

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind) {

	typedef complex<double> cd;
	typedef double d;
	if (is_same<T, cd>::value) {
		if (active1 == active2) {

			int a = active1;
			int b = before;
			int c = behind;

			rhomat_((double *) &A[0], (double *) &B[0], (double *) &S[0],
				&a, &b, &c);
			return;
		}
	}

	contraction1(S, A, B, before, active1, active2, behind, true);
}

template<typename T>
Matrix<T> contraction(const Tensor<T>& A, const Tensor<T>& B, size_t k) {
	const TensorShape& tdim_a(A.shape());
	const TensorShape& tdim_b(B.shape());
	assert(k < tdim_a.order());
	assert(k < tdim_b.order());
	size_t active1 = tdim_a[k];
	size_t active2 = tdim_b[k];
	Matrix<T> S(active1, active2);
	contraction(S, A, B, k);
	return S;
}

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k, bool zero) {
	const TensorShape& tdim_a(A.shape());
	const TensorShape& tdim_b(B.shape());
	assert(k < tdim_a.order());
	assert(k < tdim_b.order());
	size_t before = tdim_a.before(k);
	size_t after = tdim_a.after(k);
	assert(tdim_b.before(k) == before);
	assert(tdim_b.after(k) == after);
	size_t active1 = tdim_a[k];
	size_t active2 = tdim_b[k];
	assert(tdim_a.totalDimension() / active1 == tdim_b.totalDimension() / active2);
	if (zero) { S.zero(); }
	contraction(S, A, B, before, active1, active2, after);
}

template<typename T, typename U>
void matrixTensor1(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero) {

	if (zero) { C.zero(); }

	size_t dimafter = active * before;
	size_t dimafterC = activeC * before;
	for (size_t aft = 0; aft < after; ++aft) {
		for (size_t act = 0; act < active; ++act) {
			for (size_t actC = 0; actC < activeC; ++actC) {
				for (size_t bef = 0; bef < before; ++bef) {
					C[dimafterC * aft + actC * before + bef] +=
						h[act * activeC + actC] * B[dimafter * aft + act * before + bef];
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"
	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;

	if (activeB == activeC) {
		if (is_same<U, cd>::value && is_same<T, cd>::value) {
			matvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
				&a, &b, &c, &add);
			return;
		} else if (is_same<U, d>::value && is_same<T, d>::value) {
			rmatvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
				&a, &b, &c, &add);
			return;
		}
	}

	matrixTensor1(C, A, B, before, activeC, activeB, after, zero);
}

template<typename T, typename U>
void tMatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"

	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;
	if (is_same<U, cd>::value && is_same<T, cd>::value) {
		ctmatvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
			&a, &b, &c, &add);
//	} else if (is_same<U, d>::value && is_same<T, d>::value) {
//		rmatvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
//			&a, &b, &c, &add);
	} else {
		if (zero) { C.zero(); }
		size_t actbefB = activeB * before;
		size_t actbefC = activeC * before;
		size_t Cidx = 0;
		size_t Bidx = 0;
		size_t Aidx = 0;
		size_t kpreidxB = 0;
		size_t kpreidxC = 0;
		size_t Bpreidx = 0;
		size_t Cpreidx = 0;

		if (before == 1) {
//#pragma omp parallel for
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					Bidx = l + kpreidxB;
					for (size_t j = 0; j < activeC; ++j) {
						Cidx = j + kpreidxC;
						Aidx = j * activeB + l;
						C[Cidx] += conj(A[Aidx]) * B[Bidx];
//					C[Cidx] += conj(A[l, j]) * B[Bidx];
					}
				}
			}
		} else {
//#pragma omp parallel for
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					Bpreidx = l * before + kpreidxB;
					for (size_t j = 0; j < activeC; ++j) {
						Aidx = j * activeB + l;
						Cpreidx = j * before + kpreidxC;
						for (size_t i = 0; i < before; ++i) {
							Cidx = Cpreidx + i;
							Bidx = Bpreidx + i;
							C[Cidx] += conj(A[Aidx]) * B[Bidx];
//						C[Cidx] += conj(A[l, j]) * B[Bidx];
						}
					}
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape tdim(B.shape());
	TensorShape tdimC(C.shape());

	size_t after = tdim.after(mode);
	size_t before = tdim.before(mode);
	size_t active1 = A.dim1();
	size_t active2 = A.dim2();

	assert(mode < tdim.order());
	assert(A.dim2() == tdim[mode]);
	assert(A.dim1() == tdimC[mode]);

	matrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> matrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());

	if (A.dim1() == A.dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		matrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		TensorShape tdim(B.shape());
		size_t active1 = A.dim1();
		size_t active2 = A.dim2();
		tdim = replaceDimension(tdim, mode, active1);
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		assert(active1 == C.shape()[mode]);
		assert(active2 == B.shape()[mode]);
		cout << "non-quadratic mattensor implemented but tested only once so far.\n";
		matrixTensor(C, A, B, before, active1, active2, after, false);
		cout << "done.\n";
		return C;
	}
}

template<typename T, typename U>
void tensorMatrix(Tensor<T>& C, const Tensor<T>& B, const Matrix<U>& A, size_t mode, bool zero) {
	tensorMatrix(C, B, A.transpose(), mode, zero);
}

template<typename T, typename U>
Tensor<T> tensorMatrix(const Tensor<T>& B, const Matrix<U>& A, size_t mode) {
	return matrixTensor(A.transpose(), B, mode);
}

template<typename T, typename U>
Tensor<T> tMatrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	/// @TODO: remove this function; replace by tensorMatrix
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());
	assert(mode >= 0);
	assert(A.dim1() == B.shape()[mode]);

	if (A.dim1() == A.dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		tMatrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		size_t activeC = A.dim2();
		size_t activeB = A.dim1();
		TensorShape tdim(B.shape());
		tdim = replaceDimension(tdim, mode, A.dim2());
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		Tensor<T> C(tdim);
		cout << "non-quadratic mattensor implemented but not tested, yet.\n";
		tMatrixTensor(C, A, B, before, activeC, activeB, after, false);
		getchar();
		return C;
	}
}

template<typename T, typename U>
void multStateAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, bool zero) {
	const TensorShape& tdimB(B.shape());
	const TensorShape& tdimC(C.shape());

	const size_t before = tdimB.lastBefore();
	const size_t active1 = tdimB.lastDimension();
	const size_t active2 = tdimC.lastDimension();
	const size_t after = 1;

	assert(A.dim2() == active1);
	assert(A.dim1() == active2);
	assert(before == tdimC.lastBefore());

	matrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim_b(B.shape());
	size_t ntensor = tdim_b.lastDimension();
	assert(A.dim2() == ntensor);

	TensorShape tdim_c(tdim_b);
	tdim_c.setDimension(A.dim1(), tdim_c.lastIdx());
	Tensor<T> C(tdim_c);
	multStateAB(C, A, B);
	return C;
}

template<typename T, typename U>
void multStateArTB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim(B.shape());
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();
	for (size_t n = 0; n < ntensor; n++) {
		size_t B_idx = n * dimpart;
		for (size_t m = 0; m < ntensor; m++) {
			size_t C_idx = m * dimpart;
			size_t A_idx = m * ntensor;
			for (size_t i = 0; i < dimpart; i++) {
				/// C(i, m) += A(n, m) * B(i, n);
				C[C_idx + i] += A[A_idx + n] * B[B_idx + i];
			}
		}
	}
}

template<typename T, typename U>
Tensor<T> multStateArTB(const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim(B.shape());
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();
	assert(A.dim1() == A.dim2());
	assert(A.dim2() == ntensor);

	Tensor<T> C(tdim);
	multStateArTB(C, A, B);

	return C;
}

template<typename T, typename U>
void multAdd(Tensor<T>& A, const Tensor<T>& B, U coeff) {
	const TensorShape& tdim = A.shape();
	const TensorShape& tdim_2 = A.shape();
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

	TensorShape tdim(A.shape());
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
Tensor<T> qr(const Tensor<T>& A) {
	auto Amat = toMatrix(A);
	auto Qmat = qr(Amat);
	auto Q = toTensor(Qmat);
	Q.reshape(A.shape());
	return Q;
}

template<typename T>
Tensor<T> qr(const Tensor<T>& A, size_t mode) {
	auto Amat = toMatrix(A, mode);
	auto Qmat = qr(Amat);
	auto Q = toTensor(Qmat, A.shape(), mode);
	return Q;
}

//Projects B on A
template<typename T>
Tensor<T> project(const Tensor<T>& A,
	const Tensor<T>& B) {
	//calculates the overlap of A with it self
	Tensor<T> Aperp(A);
//	GramSchmidt(Aperp);
	const Matrix<T> overlap = Aperp.dotProduct(Aperp);

	//invert the overlap
	const Matrix<T> inverse_operlap = overlap.cInv();

	//calculate the scalar product of A and B
	const Matrix<T> dotproduct = Aperp.dotProduct(B);

	//multiply the scalar product and the inverse_operlap
	const Matrix<T> product = inverse_operlap * dotproduct;

	return multStateArTB(product, Aperp);
//	return multStateArTB(dotproduct, Aperp);
}

/*! \brief Project B out of A, i.e. Anew = (1-P_B)*A
 *
 * This routine takes a Tensor A and makes it orthogonal to B.
 * It can be written as A_new = (1 - P_B) A, where P_B is the
 * projector onto B.
 */
template<typename T>
Tensor<T> projectOut(const Tensor<T>& A,
	const Tensor<T>& B) {
	Tensor<T> projector = project(B, A);
	Tensor<T> perp_A(A);
	const TensorShape& tdim = A.shape();
	for (size_t i = 0; i < tdim.totalDimension(); ++i) {
		perp_A(i) -= projector(i);
	}
	return perp_A;
}

//Projects B on A
template<typename T>
Tensor<complex<double> > projectOrthogonal(const Tensor<complex<double> >& A,
	const Tensor<T>& B) {
	// calculate the scalar product of A and B
	const Matrix<complex<double> > dotproduct = A.dotProduct(B);

	return multStateArTB(dotproduct, A);
}

template<typename T>
Tensor<T> conj(Tensor<T> A) {
	for (size_t i = 0; i < A.shape().totalDimension(); ++i) {
		A[i] = conj(A[i]);
	}
	return A;
}

template<typename T>
double residual(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	auto S = A.dotProduct(A);
	return S.frobeniusNorm();
}

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A) {
	const TensorShape& shape = A.shape();
	size_t diml = shape.lastBefore();
	size_t dimr = shape.lastDimension();
	Matrix<T> B(diml, dimr);
	for (size_t j = 0; j < dimr; ++j) {
		for (size_t i = 0; i < diml; ++i) {
			B(i, j) = A(i, j);
		}
	}
	return B;
}

template<typename T>
Matrix<T> moveToMatrix(Tensor<T>& A) {
	const TensorShape& shape = A.shape();
	size_t diml = shape.lastBefore();
	size_t dimr = shape.lastDimension();
	A.ownership_ = false;
	Matrix<T> B(diml, dimr, A.coeffs_, true, false);
	return B;
}

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A, size_t mode) {
	const TensorShape& shape = A.shape();
	size_t dimbef = shape.before(mode);
	size_t dimaft = shape.after(mode);
	size_t diml = dimbef * dimaft;
	size_t dimr = shape[mode];
	Matrix<T> B(diml, dimr);
	TensorShape tmp({dimbef, dimr, dimaft});
	for (size_t bef = 0; bef < dimbef; ++bef) {
		for (size_t a = 0; a < dimr; ++a) {
			for (size_t aft = 0; aft < dimaft; ++aft) {
				size_t idxl = bef + aft * dimbef;
				size_t I = indexMapping({bef, a, aft}, tmp);
				B(idxl, a) = A(I);
			}
		}
	}
	return B;
}

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B, const TensorShape& shape, size_t mode) {
	Tensor<T> A(shape);

	size_t dimbef = shape.before(mode);
	size_t dimaft = shape.after(mode);
	size_t dimr = shape[mode];
	TensorShape tmp({dimbef, dimr, dimaft});
	for (size_t bef = 0; bef < dimbef; ++bef) {
		for (size_t a = 0; a < dimr; ++a) {
			for (size_t aft = 0; aft < dimaft; ++aft) {
				size_t idxl = bef + aft * dimbef;
				size_t I = indexMapping({bef, a, aft}, tmp);
				A(I) = B(idxl, a);
			}
		}
	}
	return A;
}

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B) {
	TensorShape shape({B.dim1(), B.dim2()});
	Tensor<T> A(shape);
	for (size_t j = 0; j < B.dim2(); ++j) {
		for (size_t i = 0; i < B.dim1(); ++i) {
			A(i, j) = B(i, j);
		}
	}
	return A;
}

template<typename T>
Tensor<T> moveToTensor(Matrix<T>& B) {
	TensorShape shape({B.dim1(), B.dim2()});
	Tensor<T> A(shape, &B[0], true, false);
	return A;
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
	if (A.shape() != B.shape()) { return false; }
	for (size_t k = 0; k < A.shape().totalDimension(); ++k) {
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
