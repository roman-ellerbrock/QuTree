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
	if (InitZero) { Zero(); }
}

template<typename T>
Tensor<T>::Tensor(const TensorShape& dim, const bool InitZero)
	:shape_(dim), coeffs_((T*)malloc(dim.totalDimension()*sizeof(T)) ), ownership_(true) {
//	:shape_(dim), coeffs_(new T[dim.totalDimension()]), ownership_(true) {
	if (InitZero) { Zero(); }
}

template<typename T>
Tensor<T>::Tensor(istream& is)
	:Tensor() {
	Read(is);
}

template<typename T>
Tensor<T>::Tensor(const string& filename)
	: Tensor() {
	ifstream is(filename);
	Read(is);
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

/*
//////////////////////////////////////////////////////////
template<typename T>
inline T& Tensor<T>::operator()(const size_t i, const size_t j, const size_t k, const size_t f, const size_t n) {
	size_t a = shape_.before(f);
	size_t b = shape_[f];
	size_t dimpart = shape_.lastBefore();
	size_t idx = n * dimpart + k * a * b + j * a + i;
	assert(i < a);
	assert(j < b);
	assert(n < shape_.lastDimension());
	assert(f < shape_.order());
	return coeffs_[idx];
}

template<typename T>
inline T& Tensor<T>::operator()(const size_t i, const size_t j, const size_t k, const size_t f, const size_t n) const {
	size_t a = shape_.before(f);
	size_t b = shape_[f];
	size_t dimpart = shape_.lastBefore();
	size_t idx = n * dimpart + k * a * b + j * a + i;
	assert(i < a);
	assert(j < b);
	assert(n < shape_.lastDimension());
	assert(f < shape_.order());
	return coeffs_[idx];
}

// Double Hole Operator
template<typename T>
T& Tensor<T>::operator()(const size_t bef, const size_t i, const size_t mid,
	const size_t j, const size_t beh, const size_t mode1, const size_t mode2,
	const size_t n) const {
	assert(mode1 < mode2);
	assert(mode2 < shape_.order());
	size_t dimpart = shape_.lastBefore();
	size_t before = shape_.before(mode1);
	size_t active1 = shape_[mode1];
	size_t active2 = shape_[mode2];
	size_t before2 = shape_.before(mode2);
	size_t idx = n * dimpart + before2 * active2 * beh + before2 * j +
		mid * before * active1 + before * i + bef;
	// mostly checking
	assert(bef < before);
	assert(i < active1);
	assert(j < active2);
	assert(n < shape_.lastDimension());
	return coeffs_[idx];
}*/

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
void Tensor<T>::Write(ostream& os) const {
	// Verification
	os.write("TENS", 4);

	// Write the TensorDim
	shape_.Write(os);

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
void Tensor<T>::Write(const string& file) const {
	ofstream os(file);
	Write(os);
}

template<typename T>
void Tensor<T>::Read(istream& is) {
	// Check if binary string contains a Tensor
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TENS");
	assert(s_key == s_check);

	// Read the TensorDim
	TensorShape newtdim;
	newtdim.ReadDim(is);

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
void Tensor<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

//////////////////////////////////////////////////////////
// Math Operators
//////////////////////////////////////////////////////////

template<typename T>
Tensor<T>& Tensor<T>::operator+=(const Tensor& A) {
	assert(A.shape().totalDimension() == shape().totalDimension());
	T const * Ax = A.coeffs_;
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
Tensor<T> Tensor<T>::AdjustDimensions(const TensorShape& newTDim) const {
	// Increase the dimensions of the Tensor from old TensorDim
	// to new TensorDim 

	assert(newTDim.order() == shape_.order());
	// Increase the active_ modes
	Tensor<T> Acoeff(*this);
	for (size_t k = 0; k < shape_.order(); k++) {
		size_t act = newTDim[k];
		Acoeff = Acoeff.AdjustActiveDim(act, k);
	}

	// Increase the number of Tensors
	size_t ntens = newTDim.lastDimension();
	Acoeff = Acoeff.AdjustStateDim(ntens);

	return Acoeff;
}

template<typename T>
Tensor<T> Tensor<T>::AdjustActiveDim(size_t active, size_t mode) const {
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
Tensor<T> Tensor<T>::AdjustStateDim(size_t n) const {
	return AdjustActiveDim(n, shape().lastIdx());
}

template<typename T>
void Tensor<T>::Reshape(const TensorShape& new_dim) {
	/// Check that total size is the same
	assert(shape_.totalDimension() == new_dim.totalDimension());
	shape_ = new_dim;
}

//////////////////////////////////////////////////////////
// Operations on Tensors
//////////////////////////////////////////////////////////
template<typename T>
T Tensor<T>::singleDotProduct(const Tensor& A, size_t n, size_t m) const {
	T result = 0;
#pragma omp parallel for reduction(+: result)
	for (size_t i = 0; i < A.shape().lastBefore(); i++) {
		result += conjugate(operator()(i, n)) * A(i, m);
	}
	return result;
}

template<typename T>
Matrix<T> Tensor<T>::DotProduct(const Tensor<T>& A) const {
	TensorShape tdima(A.shape());
	// Every tensor can have different amount of states but same dimpart

	size_t nmax = tdima.lastDimension();
	size_t mmax = shape_.lastDimension();
	size_t npart = shape_.lastBefore();
	assert(tdima.lastBefore() == npart);

	Matrix<T> S(mmax, nmax);
	Contraction(S, *this, A, shape_.lastIdx());
/*#pragma omp parallel for
	for (size_t n = 0; n < nmax; n++) {
		for (size_t m = 0; m < mmax; m++) {
			for (size_t i = 0; i < npart; i++) {
//				S(m, n) += conj(operator()(i, m))*A(i, n);
//				S(m, n) += conj(operator[](m * npart + i)) * A[n * npart + i];
				S(m, n) += conj(coeffs_[m * npart + i]) * A[n * npart + i];
			}
		}
	}
	*/
	return S;
}

template<typename T>
void Tensor<T>::Zero() {
	memset(coeffs_, 0, shape_.totalDimension() * sizeof(T));
//	for (size_t i = 0; i < shape_.totalDimension(); i++)
//		coeffs_[i] = 0;
}

template<typename T>
Tensor<T>::Tensor(const Matrix<T>& mat)
	: Tensor<T>({mat.Dim1(), mat.Dim2()}) {

	for (size_t i = 0; i < mat.Dim2(); ++i) {
		for (size_t k = 0; k < mat.Dim1(); ++k) {
			this->operator[](indexMapping({k, i}, shape_)) = mat(k, i);
		}
	}
}

// Non-member functions
//////////////////////////////////////////////////////////
template<typename T>
T SingleDotProd(const Tensor<T>& A, const Tensor<T>& B, size_t n, size_t m) {
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
void matvec_(double* C, double* B, double* mat,
	int* a, int* b, int* c, int* add);
void ctmatvec_(double* C, double* B, double* mat,
	int* a, int* b, int* c, int* add);
void rmatvec_(double* C, double* B, double* mat,
	int* a, int* b, int* c, int* add);
void rhomat_(double* Bra, double* Ket, double* M,
	int* a, int* b, int* c);
}

template<typename T>
void TensorContraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind) {

	int a = active1;
	int b = before;
	int c = behind;

	rhomat_((double*) &A[0], (double*) &B[0], (double*) &S[0],
		&a, &b, &c);

	/*
	// Variables for precalculation of indices
	size_t actbef1 = active1 * before;
	size_t actbef2 = active2 * before;
	size_t Sidx = 0;
	size_t Aidx = 0;
	size_t Bidx = 0;
	size_t ipreidx = 0;
	size_t jpreidx = 0;
	size_t npreidx1 = 0;
	size_t npreidx2 = 0;
	 */

	// Avoid unnecessary thread launches
	/*
	const char* threads = getenv("OMP_NUM_THREADS");
	if (threads) {
	    if (after < atoi(threads)) {
            omp_set_num_threads(after);
	    }
	}
	*/
	/*
#pragma omp parallel for private(npreidx1, npreidx2, jpreidx, Sidx, ipreidx, Aidx, Bidx)
	for (size_t n = 0; n < behind; n++) {
		npreidx1 = n * actbef1;
		npreidx2 = n * actbef2;
		for (size_t j = 0; j < active2; j++) {
			jpreidx = npreidx2 + j * before;
			for (size_t i = 0; i < active1; i++) {
				// S(i, j)
				Sidx = j * active1 + i;
				ipreidx = npreidx1 + i * before;
				for (size_t l = 0; l < before; l++) {
					// A(l, i, n)
					Aidx = ipreidx + l;
					// B(l, j, n)
					Bidx = jpreidx + l;
					S[Sidx] += conj(A[Aidx]) * B[Bidx];
				}
			}
		}
	}
	*/
}

template<typename T>
Matrix<T> Contraction(const Tensor<T>& A, const Tensor<T>& B, size_t k) {
	const TensorShape& tdim_a(A.shape());
	const TensorShape& tdim_b(B.shape());
	assert(k < tdim_a.order());
	assert(k < tdim_b.order());
	size_t active1 = tdim_a[k];
	size_t active2 = tdim_b[k];
	Matrix<T> S(active1, active2);
	Contraction(S, A, B, k);
	return S;
}

template<typename T>
void Contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k, bool zero) {
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
	if (zero) { S.Zero(); }
	TensorContraction(S, A, B, before, active1, active2, after);
}

template<typename T, typename U>
void MatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"
	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;

	if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
		matvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
		&a, &b, &c, &add);
	} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
		rmatvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
			&a, &b, &c, &add);
	} else {
		if (zero) { C.Zero(); }

		// Variables to Precompute index values
		size_t actbefB = activeB * before;
		size_t actbefC = activeC * before;
		size_t Aidx = 0;
		size_t Bidx = 0;
		size_t Cidx = 0;
		size_t kpreidxB = 0;
		size_t kpreidxC = 0;
		size_t lpreidx = 0;
		size_t jpreidx = 0;
		size_t lactive = 0;
		// Avoid unnecessary thread launches
		// TODO: this requires #inclue <omp.h>
		/*
		const char* threads = getenv("OMP_NUM_THREADS");
		if (threads) {
			if (after < atoi(threads)) {
				omp_set_num_threads(after);
			}
		}
		 */
		if (before == 1) {
#pragma omp parallel for private(kpreidxB, kpreidxC, Bidx, Cidx, Aidx)
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					Bidx = l + kpreidxB;
					for (size_t j = 0; j < activeC; ++j) {
						Cidx = j + kpreidxC;
						Aidx = l * activeB + j; //TODO: why is this declared twice?
						Aidx = l * activeC + j;
//					assert(Cidx < C.shape().totalDimension());
//					assert(Bidx < B.shape().totalDimension());
//					assert(Aidx < A.Dim1()*A.Dim2());
						/// C(1, j, k) += A(j, l) * B(1, l, k)
						C[Cidx] += A[Aidx] * B[Bidx];
					}
				}
			}
		} else {
#pragma omp parallel for private(Aidx, Bidx, Cidx, kpreidxB, kpreidxC, lpreidx, lactive, jpreidx)
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					lpreidx = l * before + kpreidxB;
					lactive = l * activeC;
					for (size_t j = 0; j < activeC; ++j) {
						Aidx = lactive + j;
						jpreidx = j * before + kpreidxC;
						for (size_t i = 0; i < before; ++i) {
							Cidx = jpreidx + i;
							Bidx = lpreidx + i;
//						assert(Cidx < C.shape().totalDimension());
//						assert(Bidx < B.shape().totalDimension());
//						assert(Aidx < A.Dim1()*A.Dim2());
							/// C(i, j, k) += A(j, l) * B(i, l, k)
							C[Cidx] += A[Aidx] * B[Bidx];
						}
					}
				}
			}
		}
	}
}

template<typename T, typename U>
void TMatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"

	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;
	if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
		ctmatvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
			&a, &b, &c, &add);
//	} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
//		rmatvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
//			&a, &b, &c, &add);
	} else {
		if (zero) { C.Zero(); }
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
void MatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape tdim(B.shape());
	TensorShape tdimC(C.shape());

	size_t after = tdim.after(mode);
	size_t before = tdim.before(mode);
	size_t active1 = A.Dim1();
	size_t active2 = A.Dim2();

	assert(mode < tdim.order());
	assert(A.Dim2() == tdim[mode]);
	assert(A.Dim1() == tdimC[mode]);

	MatrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> MatrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());

	if (A.Dim1() == A.Dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		MatrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		TensorShape tdim(B.shape());
		size_t active1 = A.Dim1();
		size_t active2 = A.Dim2();
		tdim = replaceDimension(tdim, mode, active1);
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		assert(active1 == C.shape()[mode]);
		assert(active2 == B.shape()[mode]);
		cout << "non-quadratic mattensor implemented but tested only once so far.\n";
		MatrixTensor(C, A, B, before, active1, active2, after, false);
		return C;
	}
}

template<typename T, typename U>
void TensorMatrix(Tensor<T>& C, const Tensor<T>& B, const Matrix<U>& A, size_t mode, bool zero) {
	TensorMatrix(C, B, A.Transpose(), mode, zero);
}

template<typename T, typename U>
Tensor<T> TensorMatrix(const Tensor<T>& B, const Matrix<U>& A, size_t mode) {
	return MatrixTensor(A.Transpose(), B, mode);
}

template<typename T, typename U>
Tensor<T> multATB(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());
	assert(mode >= 0);
	assert(A.Dim1() == B.shape()[mode]);

	if (A.Dim1() == A.Dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		TMatrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		size_t activeC = A.Dim2();
		size_t activeB = A.Dim1();
		TensorShape tdim(B.shape());
		tdim = replaceDimension(tdim, mode, A.Dim2());
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		Tensor<T> C(tdim);
		cout << "non-quadratic mattensor implemented but not tested, yet.\n";
		TMatrixTensor(C, A, B, before, activeC, activeB, after, false);
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

	assert(A.Dim2() == active1);
	assert(A.Dim1() == active2);
	assert(before == tdimC.lastBefore());

	MatrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim_b(B.shape());
	size_t ntensor = tdim_b.lastDimension();
	assert(A.Dim2() == ntensor);

	TensorShape tdim_c(tdim_b);
	tdim_c.setDimension(A.Dim1(), tdim_c.lastIdx());
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
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim2() == ntensor);

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
void GramSchmidt(Tensor<T>& A) {
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
				T overlap = SingleDotProd(A, A, m, n);
				accumoverlap += abs(overlap);
				for (size_t i = 0; i < dimpart; i++) {
					A(i, n) -= overlap * A(i, m);
				}
			}

			// renormalize
			T norm = SingleDotProd(A, A, n, n);
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

Tensorcd QR(const Tensorcd& A) {
	auto Amat = toMatrix(A);
	auto Qmat = QR(Amat);
	auto Q = toTensor(Qmat);
	Q.Reshape(A.shape());
	return Q;
}


//Projects B on A
template<typename T>
Tensor<complex<double> > Project(const Tensor<complex<double> >& A,
	const Tensor<T>& B) {
	//calculates the overlap of A with it self
	Tensor<complex<double> > Aperp(A);
//	GramSchmidt(Aperp);
	const Matrix<complex<double> > overlap = Aperp.DotProduct(Aperp);

	//invert the overlap
	const Matrix<complex<double> > inverse_operlap = overlap.cInv();

	//calculate the scalar product of A and B
	const Matrix<complex<double> > dotproduct = Aperp.DotProduct(B);

	//multiply the scalar product and the inverse_operlap
	const Matrix<complex<double> > product = inverse_operlap * dotproduct;

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
Tensor<T> ProjectOut(const Tensor<T>& A,
	const Tensor<T>& B) {
	Tensorcd projector = Project(B, A);
	Tensorcd perp_A(A);
	const TensorShape& tdim = A.shape();
	for (size_t i = 0; i < tdim.totalDimension(); ++i) {
		perp_A(i) -= projector(i);
	}
	return perp_A;
}

//Projects B on A
template<typename T>
Tensor<complex<double> > ProjectOrthogonal(const Tensor<complex<double> >& A,
	const Tensor<T>& B) {
	// calculate the scalar product of A and B
	const Matrix<complex<double> > dotproduct = A.DotProduct(B);

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
double Residual(Tensor<T> D, const Tensor<T>& B) {
	D -= B;
	auto S = D.DotProduct(D);
	return S.FrobeniusNorm();
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
Tensor<T> toTensor(const Matrix<T>& B) {
	TensorShape shape({B.Dim1(), B.Dim2()});
	Tensor<T> A(shape);
	for (size_t j = 0; j < B.Dim2(); ++j) {
		for (size_t i = 0; i < B.Dim1(); ++i) {
			A(i, j) = B(i, j);
		}
	}
	return A;
}

template<typename T>
ostream& operator<<(ostream& os, const Tensor<T>& A) {
	A.Write(os);
	return os;
}

template<typename T>
istream& operator>>(istream& is, Tensor<T>& A) {
	A.Read(is);
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
