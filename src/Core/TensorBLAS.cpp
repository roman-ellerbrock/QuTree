//
// Created by Roman Ellerbrock on 5/22/21.
//

#include "Core/TensorBLAS.h"
#include <cblas.h>
#include <numeric>
#include "Core/MatrixBLAS.h"


#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta) {
	// A[dim1, dim2] --> A[dim2, dim1]
	/// simple in-place transpose
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
//			dest[j + dim2 * i] = src[i + dim1 * j];
			dest[j + dim2 * i] = beta * dest[j + dim2 * i] + src[i + dim1 * j];
		}
	}
}

template<typename T, int blocksize>
void transpose2(T *dest, const T *src, size_t lda, size_t ldb) {
	/// Incorporate blocking to minimize cache misses
	// dest[b, a, c] = src[a, b, c]
	/// simple in-place transpose
	size_t nblockA = lda / blocksize;
	size_t nblockB = ldb / blocksize;
	for (size_t b = 0; b < nblockB; ++b) {
		for (size_t a = 0; a < nblockA; ++a) {
			for (size_t bb = 0; bb < blocksize; ++bb) {
				for (size_t aa = 0; aa < blocksize; ++aa) {
					// dest[a', b'] <-- src[b', a']
					// a' = a * blocksize + aa
					// b' = b * blocksize + bb
					// idxDest = b' * lda + a' = (b * blocksize + bb) * lda + a * blocksize + aa
					// idxSrc  = a' * ldb + b' = (a * blocksize + aa) * ldb + b * blocksize + bb
					size_t idxDest = (b * blocksize + bb) * lda + a * blocksize + aa;
					size_t idxSrc = (a * blocksize + aa) * ldb + b * blocksize + bb;
					dest[idxDest] = src[idxSrc];
				}
			}
		}
	}

	// right block and square block
	size_t restA = lda % blocksize;
	size_t offsetA = lda - restA;
	for (size_t a = 0; a < restA; ++a) {
		for (size_t b = 0; b < ldb; ++b) {
			// dest[a', b'] <-- src[a', b']
			size_t aa = offsetA + a;
			dest[aa + b * lda] = src[aa * ldb + b];
		}
	}
	// bottom block
	size_t restB = lda % blocksize;
	size_t offsetB = ldb - restB;
	for (size_t b = 0; b < restB; ++b) {
		for (size_t a = 0; a < nblockA * blocksize; ++a) {
			// dest[a', b'] <-- src[a', b']
			size_t bb = offsetB + b;
			dest[a + bb * lda] = src[a * ldb + bb];
		}
	}
}

template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[b, a, c] = src[a, b, c]
	for (size_t c = 0; c < C; ++c) {
		//	transpose2<T,4>(&dest[c * A * B], &src[c * A * B], A, B);
		transpose(&dest[c * A * B], &src[c * A * B], A, B);
	}
}

template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[a, c, b] = src[a, b, c]
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			memcpy(&dest[c * A + b * A * C], &src[b * A + c * A * B], A * sizeof(T));
/*			for (size_t a = 0; a < A; ++a) {
				dest[a + c * A + b * A * C] = src[a + b * A + c * A * B];
			}*/
		}
	}
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
void matrixTensor2(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	Tensor<T>& D, size_t before, size_t active, size_t activeC, size_t after, bool zero) {
	/// D is work tensor with shape of C.
	typedef complex<double> cd;
	typedef double d;
	/// If not double/complex<double> type, fall back to straightforward implementation
	if constexpr(!((is_same<U, cd>::value && is_same<T, cd>::value))
		|| (is_same<U, d>::value && is_same<T, d>::value)) {
		matrixTensor1(C, h, B, before, active, activeC, after, zero);
		return;
	}

	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }

	if (before == 1) {
		size_t m = activeC;
		size_t k = active; //activeB
		size_t n = after;
		if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
				(void *) &z, (void *) &h[0], m, (void *) &B[0], k, (void *) &zz, (void *) &C[0], m);
		} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
				z, (double *) &h[0], m, (double *) &B[0], k, zz, (double *) &C[0], m);
		}
		return;
	}

	size_t m = activeC;
	size_t k = active; //activeB
	size_t n = before;

	size_t pref = before * active;
	size_t prefC = before * activeC;
	T zer = 0.;
	for (size_t aft = 0; aft < after; ++aft) {
		if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
				(void *) &z, (void *) &h[0], m, (void *) &B[aft * pref], n, (void *) &zer, (void *) &D[aft * prefC], m);
		} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
				z, (double *) &h[0], m, (double *) &B[aft * pref], n, zer, (double *) &D[aft * prefC], m);
		}
//		transpose2<T, 4>(&C[aft * prefC], &D[aft * prefC], activeC, before);
		transpose(&C[aft * prefC], &D[aft * prefC], activeC, before, zz);
	}
	if (!zero) { }
}

template<typename T, typename U>
void matrixTensor3(Tensor<T>& hKet, const Matrix<U>& h, const Tensor<T>& Ket,
	Tensor<T>& Ket_work, Tensor<T>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	/// hKet[a, b2, c] = h[b2, b] * Ket[a, b, c]
	/// A, B, B2, C are dimensions of indices a, b, b2, c
	/// Transpose tensor and map to BLAS ?geem
	typedef complex<double> cd;
	typedef double d;

	size_t AC = A * C;
	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	transposeAB(&Ket_work[0], &Ket[0], A, B, C);
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, B2, AC, B,
		(void *) &z, (void *) &h[0], B2, (void *) &Ket_work[0], B,
		(void *) &zz, (void *) &hKet_work[0], B2);

	transposeAB(&hKet[0], &hKet_work[0], B2, A, C);
}

//====================================================================

//====================================================================

template<typename T>
void contraction2(Matrix<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	Tensor<T>& bra_work, Tensor<T>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	typedef complex<double> cd;
	typedef double d;
	/// If not double/complex<double> type, fall back to straightforward implementation
	if constexpr(!(is_same<T, cd>::value || is_same<T, d>::value)) {
		contraction1(h, bra, ket, A, B, B2, C, zero);
		return;
	}

	transposeBC(&bra_work[0], &bra[0], A, B, C);
	transposeBC(&ket_work[0], &ket[0], A, B2, C);

	size_t AC = A * C;
	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	if constexpr(is_same<T, cd>::value) {
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
			(void *) &z, (void *) &bra_work[0], AC, (void *) &ket_work[0], AC, (void *) &zz, (void *) &h[0], B);
	} else if constexpr(is_same<T, d>::value) {
		cblas_dgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
			z, (double *) &bra_work[0], AC, (double *) &ket_work[0], AC, zz, (double *) &h[0], B);
	}
}

/// ========================================================================
/// Wrappers for matrixTensor product and Tensor contraction
/// ========================================================================

/// Wrapper for matrix-Tensor product
template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, Tensor<T>& workC, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {

	TensorShape tdim(B.shape());
	TensorShape tdimC(C.shape());

	if (mode >= tdim.order()) {
		cerr << "matrixTensor error: mode too large.\n";
		exit(1);
	}

	size_t after = tdim.after(mode);
	size_t before = tdim.before(mode);
	size_t active1 = A.dim1();
	size_t active2 = A.dim2();

	if (!(A.dim2() == tdim[mode])) {
		cerr << "matrix Tensor error: active dimension wrong.\n";
		exit(1);
	}
	if (!(A.dim1() == tdimC[mode])) {
		cerr << "matrix Tensor error: left active dimension wrong.\n";
		exit(1);
	}
	if (tdim.before(mode) != tdimC.before(mode)) {
		cerr << "matrix Tensor error: before dimension wrong.\n";
		exit(1);
	}
	if (tdim.after(mode) != tdimC.after(mode)) {
		cerr << "matrix Tensor error: after dimension wrong.\n";
		exit(1);
	}
	if (tdimC.totalDimension() != workC.shape().totalDimension()) {
		cerr << "matrix Tensor error: work array has wrong dimension.\n";
		exit(1);
	}

	matrixTensor2(C, A, B, workC, before, active1, active2, after, zero);
}

template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	Tensor<T> workC(C);
	matrixTensorBLAS(C, workC, A, B, mode, zero);
}

template<typename T>
void contractionBLAS(Matrix<T>& h, Tensor<T>& workA, Tensor<T>& workB, const Tensor<T>& A, const Tensor<T>& B,
	size_t mode, bool zero) {

	TensorShape tdimA(A.shape());
	TensorShape tdimB(B.shape());
	if (mode >= tdimA.order()) {
		cerr << "contraction error: mode too large for Bra.\n";
		exit(1);
	}
	if (mode >= tdimB.order()) {
		cerr << "contraction error: mode too large for Ket.\n";
		exit(1);
	}

	size_t after = tdimA.after(mode);
	size_t before = tdimA.before(mode);
	size_t activeA = tdimA[mode];
	size_t activeB = tdimB[mode];

	if (!(h.dim1() == activeA)) {
		cerr << "contraction error: ket active dimension wrong.\n";
		exit(1);
	}
	if (!(h.dim2() == activeB)) {
		cerr << "contraction error: bra active dimension wrong.\n";
		exit(1);
	}
	if (tdimA.before(mode) != tdimB.before(mode)) {
		cerr << "contraction error: before dimension wrong.\n";
		exit(1);
	}
	if (tdimA.after(mode) != tdimB.after(mode)) {
		cerr << "contraction error: after dimension wrong.\n";
		exit(1);
	}

	if (tdimA.totalDimension() != workA.shape().totalDimension()) {
		cerr << "contraction error: work array A has wrong dimension.\n";
		exit(1);
	}
	if (tdimB.totalDimension() != workB.shape().totalDimension()) {
		cerr << "contraction error: work array B has wrong dimension.\n";
		exit(1);
	}

	contraction2(h, A, B, workA, workB, before, activeA, activeB, after, zero);
}

/// Wrapper for matrix-Tensor product
template<typename T>
void contractionBLAS(Matrix<T>& h, const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero) {
	Tensor<T> Awork(A.shape());
	Tensor<T> Bwork(B.shape());
	contractionBLAS(h, Awork, Bwork, A, B, mode, zero);
}

template<typename T>
Matrix<T> contractionBLAS(const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape tdimA(A.shape());
	TensorShape tdimB(B.shape());
	if (mode >= tdimA.order()) {
		cerr << "contraction error: mode too large for Bra.\n";
		exit(1);
	}
	if (mode >= tdimB.order()) {
		cerr << "contraction error: mode too large for Ket.\n";
		exit(1);
	}

	size_t dim1 = tdimA[mode];
	size_t dim2 = tdimB[mode];
	Matrix<T> h(dim1, dim2);
	contractionBLAS(h, A, B, mode, zero);
	return h;
}

template<typename T, typename U>
Tensor<T> matrixTensorBLAS(const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape Cshape = B.shape();
	Cshape.setDimension(A.dim1(), mode);
	Tensor<T> C(Cshape);
	matrixTensorBLAS(C, A, B, mode, zero);
	return C;
}

void dgeem(Matrixd& h, const Matrixd& bra, const Matrixd& ket) {

	double z = 1.0, zz = 0.;
	// test 1
	size_t m = h.dim1();
	size_t n = h.dim2();
	size_t k = bra.dim2();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
		z, &bra[0], m, &ket[0], n, zz, &h[0], m);
	// test 2
/*	size_t m = h.dim1();
	size_t n = h.dim2();
	size_t k = bra.dim1();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
		z, (double *) &bra[0], k, (double *) &ket[0], k, zz, (double *) &h[0], m);
		*/
}

// Tensor Contraction along arbitrary indices
// C_{a,c,b,d,e} = Sum_{f,g,h} A_{f,g,a,h,c}^{*} x B_{g,f,h,b,d,e}
// the general order of indices is conserved
//
// the input arguments are the starting tensors A,B
// the result tensor result
// and the indices along which the contraction will be done
// for the example above, the call would look like:
// general_contraction(A, B, C, {0,1,3},{1,0,4}),
// meaning index 0 of A will be contracted with index 1 of B etc.


template<typename T>
void general_contraction(const Tensor<T>& A,
                         const Tensor<T>& B,
                         Tensor<T>& result,
                         const vector<size_t> &A_indices,
                         const vector<size_t> &B_indices){

    // is contraction legal?
    is_contraction_legal(A,B,A_indices,B_indices);

    // now the contractions can be started
    // to become efficient, the tensors will be transposed until they are ordered correctly
    // e.g.: A_indices = {3,1,2}, B_indices = {3,6,1},
    // then the indices will be ordered such that the indices {3,1,2} are the first ones of A,
    // and the indices {3,6,1} are the first ones of B


    // now transpose A first
    // for this, first transpose the relevant (to-be-contracted) indices to the front,
    // then transpose the irrelevant indices into the correct order for the result
    auto A_correct_form{A};
    auto B_correct_form{B};
    general_transpose_to_order(A_correct_form,A,A_indices);
    general_transpose_to_order(B_correct_form,B,B_indices);

    // these two tensors are working objects now, calculate correct and intermediate forms
    // for the tensor contraction and the result form
    std::vector<size_t> result_form;

    // calculate size of first (super)-index along which the contraction will be done
    size_t A_contraction_length = 1;
    for(size_t i = 0; i < A_indices.size(); ++i){
        A_contraction_length *= A_correct_form.shape().dimensions()[i];
    }

    size_t B_contraction_length = 1;
    for(size_t i = 0; i < B_indices.size(); ++i){
        B_contraction_length *= B_correct_form.shape().dimensions()[i];
    }


    // set other sizes for the result tensor and the intermediate tensors
    size_t A_second_index = 1;
    for(size_t i = A_indices.size(); i < A.shape().order(); ++i){
        result_form.push_back(A_correct_form.shape().dimensions()[i]);
        A_second_index *= A_correct_form.shape().dimensions()[i];

    }
    size_t B_second_index = 1;
    for(size_t i = B_indices.size(); i < B.shape().order(); ++i){
        result_form.push_back(B_correct_form.shape().dimensions()[i]);
        B_second_index *= B_correct_form.shape().dimensions()[i];
    }

    TensorShape new_A_shape({A_contraction_length,A_second_index});
    TensorShape new_B_shape({B_contraction_length,B_second_index});
    TensorShape result_shape(result_form);

    A_correct_form.shape() = new_A_shape;
    B_correct_form.shape() = new_B_shape;

    // now the contraction can be performed
    Tensor<T> result_tensor{result_shape};

    // do the contraction itself

    if constexpr(is_same<T, std::complex<double>>::value) {
        const std::complex<double> one{1.,0.};
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
                    A_second_index, // rows of A
                    B_second_index, // cols of B
                    B_contraction_length, // cols of A = rows of B
                    (void*)&one,
                    (void*)&(A_correct_form.coeffs_[0]),
                    B_contraction_length,
                    (void*)&(B_correct_form.coeffs_[0]),
                    B_contraction_length,
                    (void*)&one,
                    (void*)&(result_tensor.coeffs_[0]),
                    A_second_index);

    } else if constexpr(is_same<T, double>::value) {
        const double one = 1.;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    A_second_index, // rows of A
                    B_second_index, // cols of B
                    B_contraction_length, // cols of A = rows of B
                    one,
                    &(A_correct_form.coeffs_[0]),
                    B_contraction_length,
                    &(B_correct_form.coeffs_[0]),
                    B_contraction_length,
                    one,
                    &(result_tensor.coeffs_[0]),
                    A_second_index);
    }

    result = std::move(result_tensor);
}

// helper function determining if a tensor requested contraction is legal
// a contraction is valid iff (i)   every index is only contracted once (present once in each vector)
//                            (ii)  the smallest index is 0
//                            (iii) the biggest index is <= rank of tensor
//                            (iv)  both contractions do have the same size
//                            (v)   do contraction dimensions match up
// all of this can be determined by sorting the indices and checking the values
template<typename T>
bool is_contraction_legal(const Tensor<T> &TensorA,
                          const Tensor<T> &TensorB,
                          vector<size_t> Acontraction,
                          vector<size_t> Bcontraction) {


    // check if the contraction index numbers match
    if(Acontraction.size() != Bcontraction.size()){
        return false;
    }

    // check if contraction dimensions match up
    for(size_t i = 0; i < Acontraction.size(); ++i){
        if(TensorA.shape().dimensions()[Acontraction[i]] != TensorB.shape().dimensions()[Bcontraction[i]]){
            return false;
        }
    }

    // sort index data
    std::sort(Acontraction.begin(), Acontraction.end());
    std::sort(Bcontraction.begin(), Bcontraction.end());


    // check whether the left and right contractions are legal on their own
    const bool A_shape_is_ok = !Acontraction.empty()
            and (Acontraction[0] >= 0)
            and (Acontraction.back() <= TensorA.shape().order());

    const bool B_shape_is_ok = !Acontraction.empty()
                         and (Acontraction[0] >= 0)
                         and (Acontraction.back() <= TensorA.shape().order());


    // check if any index is doubly present
    const bool A_doubly_present = (std::unique(Acontraction.begin(), Acontraction.end()) == Acontraction.end());
    const bool B_doubly_present = (std::unique(Bcontraction.begin(), Bcontraction.end()) == Bcontraction.end());

    // if everything is true, then the contraction is legal
    return A_shape_is_ok and B_shape_is_ok and A_doubly_present and B_doubly_present;
}

// this function transforms the src tensor to a new form specified by the input argument form
// e.g. src is a {1,2,3,4,5} tensor, when form contains {2,1,3} then the resulting tensor is
// has the shape: {3,2,4,1,5}, so the order of all other indices is preserved
// DANGER: this function does _not_ check whether form is a correct form!!
template<typename T>
void general_transpose_to_order(Tensor<T>& dst, const Tensor<T>& src, const vector<size_t> &form){

    // working tensors
    Tensor<T> tmp1{src};
    dst = src;

    // save old index order
    std::vector<int> current_order(src.shape().order());
    std::iota(current_order.begin(), current_order.end(), 0);

    // transpose the relevant indices to the correct spots (from left to right)
    int current_index = 0;
    for(const auto& i : form){

        // save last result to tmp1 again such that dst will be the result after this iteration
        if(current_index != 0){
            swap(Tensor<T>, tmp1, dst);
        }

        // first: find where this index is currently
        int current_place = -1;
        for(int j = 0; j < current_order.size(); ++j){
            if(i == current_order[j]){
                current_place = j;
                break;
            }
        }

        // now transpose properly
        general_transpose(dst,tmp1,current_index,current_place);
        swap(int, current_order[current_index], current_order[current_place]);

        current_index++;
    }

    // now reorder all other indices so they are in the correct order (ascending)
    while(current_index < src.shape().order()){
        // save last result to tmp1 again such that dst will be the result after this iteration
        swap(Tensor<T>, tmp1, dst);

        // find current smallest index and position
        int smallest = std::numeric_limits<int>::max();
        int position = 0;
        for(int i = current_index; i < src.shape().order(); ++i){
            if(current_order[i] < smallest){
                smallest = current_order[i];
                position = i;
            }
        }

        // transpose accordingly
        general_transpose(dst,tmp1,current_index,position);
        swap(int, current_order[current_index], current_order[position]);

        current_index++;
    }
}

// this is a helper function which transposes two given indices of a tensor
template<typename T>
void general_transpose(Tensor<T>& dst, const Tensor<T>& src, size_t index_one, size_t index_two){
    // prepare transpose
    dst = src;

    // check for trivial cases
    if(index_one == index_two){
        return;
    }

    // check if transpose is allowed
    assert(index_one>=0);
    assert(index_two>=0);
    assert(index_one<=src.shape().order());
    assert(index_two<=src.shape().order());

    // determine a,b,c,d,e for 5-index-transpose
    // make index_one < index_two
    if(index_one > index_two) {
        swap(size_t, index_one, index_two); // std::swap would be nicer, but macro on top of this file...
    }
    size_t b = src.shape().dimensions()[index_one];
    size_t d = src.shape().dimensions()[index_two];

    size_t a = 1;
    for(size_t i = 0; i < index_one; ++i){
        a *= src.shape().dimensions()[i];
    }
    size_t c = 1;
    for(size_t i = index_one + 1; i < index_two; ++i){
        c *= src.shape().dimensions()[i];
    }
    size_t e = 1;
    for(size_t i = index_two + 1; i < src.shape().order(); ++i){
        e *= src.shape().dimensions()[i];
    }

    general_transpose_bd(dst.coeffs_,src.coeffs_,a,b,c,d,e);

    // determine and set new shape
    auto old_dimensions = dst.shape().dimensions();
    swap(int, old_dimensions[index_one], old_dimensions[index_two]);
    TensorShape new_shape(old_dimensions);
    dst.shape() = new_shape;

}


#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
// this is a helper-function doing a 5-index-transpose
template<typename T>
void general_transpose_bd(T* dst, const T* src, size_t a, size_t b, size_t c, size_t d, size_t e) {
    // src(a,b,c,d,e) -> dst(a,d,c,b,e)
    // primitive algorithm, smarter is definitely possible

    // todo more(?) or better specializations

    // first special case: c == 1, then it becomes multiple applications of the 3-index-transpose
    if (c == 1 and a != 1) { // this is an okay implementation, even for e != 1
        const size_t d_e = a * d * b;
        const size_t s_e = a * b * d;
        for (size_t i_e = 0; i_e < e; ++i_e) {
            transposeBC(&dst[i_e * d_e], &src[i_e * s_e], a, b, d);
        }
    } else if(c == 1 and a == 1) { // second special case: if c = a = 1, then it is just another 3-index transpose
        const size_t d_e = d * b;
        const size_t s_e = b * d;
        transposeAB(dst, src, b, d, e);
    } else { // general case: slow but honest work

        // strides of src
        const size_t s_a = 1;
        const size_t s_b = a;
        const size_t s_c = a * b;
        const size_t s_d = a * b * c;
        const size_t s_e = a * b * c * d;

        // strides of dst
        const size_t d_a = 1;
        const size_t d_b = a;
        const size_t d_c = a * d;
        const size_t d_d = a * d * c;
        const size_t d_e = a * d * c * b;

        // slow but honest work
//#pragma omp parallel for collapse(5)
        for(size_t i_e = 0; i_e < e; ++i_e){
            for(size_t i_d = 0; i_d < d; ++i_d){
                for(size_t i_c = 0; i_c < c; ++i_c){
                    for(size_t i_b = 0; i_b < b; ++i_b){
                        for(size_t i_a = 0; i_a < a; ++i_a){

                            const size_t src_idx = i_a * s_a
                                                   + i_b * s_b
                                                   + i_c * s_c
                                                   + i_d * s_d
                                                   + i_e * s_e;

                            const size_t dst_idx = i_a * d_a
                                                   + i_d * d_b
                                                   + i_c * d_c
                                                   + i_b * d_d
                                                   + i_e * d_e;

                            dst[dst_idx] = src[src_idx];
                        }
                    }
                }
            }
        }

    }
}
#pragma clang diagnostic pop

// see documentation of other general contraction
// this function takes instead of A_indices and B_indices
// pairs of indices to be contracted
// it is just a wrapper
template<typename T, typename Q>
void general_contraction(const Tensor<T>& A,
                         const Tensor<T>& B,
                         Tensor<T>& result,
                         const std::vector<std::pair<Q,Q>>& contraction_pairs){

    std::vector<size_t> A_indices;
    std::vector<size_t> B_indices;

    for(const auto i : contraction_pairs){
        A_indices.push_back(i.first);
        B_indices.push_back(i.second);
    }
    general_contraction(A,B,result,A_indices,B_indices);
}




/// ========================================================================
/// Template Instantiations
/// ========================================================================

typedef complex<double> cd;

typedef double d;

template void matrixTensor1(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor2(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	Tensorcd& D, size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor3(Tensor<cd>& hKet, const Matrix<cd>& h, const Tensor<cd>& Ket,
	Tensor<cd>& Ket_work, Tensor<cd>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void contraction2(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	Tensor<cd>& bra_work, Tensor<cd>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2, cd beta);

template void transpose2<complex<double>, 4>(cd *dest, const cd *src, size_t lda, size_t ldb);

template void transposeAB(cd *dest, const cd *src, size_t A, size_t B, size_t C);

template void contraction2(Matrix<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	Tensor<d>& bra_work, Tensor<d>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void general_contraction(const Tensor<cd>& A, const Tensor<cd>& B, Tensor<cd>& result,
                                  const vector<size_t> &A_indices, const vector<size_t> &B_indices);

template void general_contraction(const Tensor<d>& A, const Tensor<d>& B, Tensor<d>& result,
                                  const vector<size_t> &A_indices, const vector<size_t> &B_indices);

template void general_transpose_bd(cd* dst, const cd* src, size_t a, size_t b, size_t c, size_t d, size_t e);

template void general_transpose_bd(d* dst, const d* src, size_t a, size_t b, size_t c, size_t d, size_t e);

template bool is_contraction_legal(const Tensor<cd> &TensorA,
                                   const Tensor<cd> &TensorB,
                                   vector<size_t> Acontraction,
                                   vector<size_t> Bcontraction);

template bool is_contraction_legal(const Tensor<d> &TensorA,
                                   const Tensor<d> &TensorB,
                                   vector<size_t> Acontraction,
                                   vector<size_t> Bcontraction);

template void general_transpose(Tensor<cd >& dst, const Tensor<cd>& src, size_t index_one, size_t index_two);
template void general_transpose(Tensor<d>& dst, const Tensor<d>& src, size_t index_one, size_t index_two);

template void general_contraction(const Tensor<cd>& A,
                         const Tensor<cd>& B,
                         Tensor<cd>& result,
                         const std::vector<std::pair<int,int>>& contraction_pairs);

template void general_contraction(const Tensor<cd>& A,
                                  const Tensor<cd>& B,
                                  Tensor<cd>& result,
                                  const std::vector<std::pair<size_t,size_t>>& contraction_pairs);

template void general_contraction(const Tensor<d>& A,
                                  const Tensor<d>& B,
                                  Tensor<d>& result,
                                  const std::vector<std::pair<int,int>>& contraction_pairs);

template void general_contraction(const Tensor<d>& A,
                                  const Tensor<d>& B,
                                  Tensor<d>& result,
                                  const std::vector<std::pair<size_t,size_t>>& contraction_pairs);

/// ==== Wrappers ====
// complex double
template void matrixTensorBLAS(Tensor<cd>& C, Tensor<cd>& workC, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void matrixTensorBLAS(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<cd>& h, Tensor<cd>& workA, Tensor<cd>& workB, const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<cd>& h, const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> matrixTensorBLAS(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Matrix<cd> contractionBLAS(const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
// double
template void matrixTensorBLAS(Tensor<d>& C, Tensor<d>& workC, const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void matrixTensorBLAS(Tensor<d>& C, const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<d>& h, Tensor<d>& workA, Tensor<d>& workB, const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<d>& h, const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template Tensor<d> matrixTensorBLAS(const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template Matrix<d> contractionBLAS(const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
