//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Tensor/TensorBLAS2.h"
#include <errno.h>

#define CHECK(X) ({int __val = (X); __val == -1 ? \
    ({ fprintf(stderr, "ERROR (" __FILE__ ":%d) -- %s\n",__LINE__,strerror(errno)); \
    exit(-1);-1;}) : __val; })


typedef complex<double> cd;
typedef double d;

size_t nrows(const TensorShape& shape, blas::Op op) {
	size_t n = shape.lastBefore();
	if (op != blas::Op::NoTrans) { n = shape.lastDimension(); }
	return n;
}

size_t ncols(const TensorShape& shape, blas::Op op) {
	size_t n = shape.lastDimension();
	if (op != blas::Op::NoTrans) { n = shape.lastBefore(); }
	return n;
}

template<typename T>
void gemm(Tensor<T>& c, const Tensor<T>& a, const Tensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b) {

	/// m: Number of rows of the matrix C and \(op(A)\). m >= 0.
	size_t m = nrows(a.shape_, op_a);
	/// @TODO: check nrows(c) == m

	/// n: Number of cols of the matrix C and \(op(B)\). n >= 0.
	size_t n = ncols(b.shape_, op_b);
	/// @TODO: check ncols(c) == n

	size_t k = ncols(a.shape_, op_a);
	/// @TODO: check ncols(c) == k

	size_t lda = (op_a == blas::Op::NoTrans) ? m : k;
	size_t ldb = (op_b == blas::Op::NoTrans) ? k : n;
	size_t ldc = m;

	blas::gemm(blas::Layout::ColMajor, op_a, op_b, m, n, k,
		alpha, a.coeffs_, lda, b.coeffs_, ldb, beta, c.coeffs_, ldc);
}

template void gemm(Tensor<cd>& c, const Tensor<cd>& a, const Tensor<cd>& b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b);
template void gemm(Tensor<d>& c, const Tensor<d>& a, const Tensor<d>& b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b);


template<typename T>
Tensor<T> gemm(const Tensor<T>& a, const Tensor<T>& b,
	T alpha, blas::Op op_a, blas::Op op_b) {
	size_t m = nrows(a.shape_, op_a);
	size_t n = ncols(b.shape_, op_b);

	Tensor<T> c({m, n});
	T beta = 0.;
	gemm(c, a, b, alpha, beta, op_a, op_b);
	return c;
}

template Tensor<cd> gemm(const Tensor<cd>& a, const Tensor<cd>& b,
	cd alpha, blas::Op op_a, blas::Op op_b);
template Tensor<d> gemm(const Tensor<d>& a, const Tensor<d>& b,
	d alpha, blas::Op op_a, blas::Op op_b);

template<typename T>
void gemmRef(Tensor<T>& c, Tensor<T> a, Tensor<T> b,
	T alpha, T beta, blas::Op op_a, blas::Op op_b) {

	if (op_a == blas::Op::Trans) {
		a = transpose(a);
	} else if(op_a == blas::Op::ConjTrans) {
		a = adjoint(a);
	}

	if (op_b == blas::Op::Trans) {
		b = transpose(b);
	} else if(op_b == blas::Op::ConjTrans) {
		b = adjoint(b);
	}

	size_t m = a.shape_[0];
	size_t k = a.shape_[1];
	size_t n = b.shape_[1];
	for (size_t M = 0; M < m; ++M) {
		for (size_t N = 0; N < n; ++N) {
			c(M, N) = beta * c(M, N);
			for (size_t K = 0; K < k; ++K) {
				c(M, N) += alpha * a(M, K) * b(K, N);
			}
		}
	}

}

template void gemmRef(Tensor<cd>& c, Tensor<cd> a, Tensor<cd> b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b);

template void gemmRef(Tensor<d>& c, Tensor<d> a, Tensor<d> b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b);

template<typename T, typename U>
void matrixTensorRef(Tensor<T>& C, const Tensor<U>& h, const Tensor<T>& B,
	size_t k, bool zero) {
	/// This is a simple implementation of the matrix-tensor product
	/// which is not optimized for efficiency.
	const TensorShape& shape = B.shape_;
	size_t before = shape.before(k);
	size_t active = shape[k];
	size_t activeC = C.shape_[k];
	size_t after = shape.after(k);
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

template void matrixTensorRef(Tensor<d>& C, const Tensor<d>& h, const Tensor<d>& B,
	size_t k, bool zero);
template void matrixTensorRef(Tensor<cd>& C, const Tensor<cd>& h, const Tensor<cd>& B,
	size_t k, bool zero);

template<typename T>
void matrixTensor(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	size_t k1, T alpha, T beta, blas::Op op_h) {
	/// D is work tensor with shape of C.
	const TensorShape& shape = B.shape_;
	size_t before = shape.before(k1);
	size_t active = shape[k1];
	size_t activeC = C.shape_[k1];
	size_t after = shape.after(k1);

	blas::Layout layout = blas::Layout::ColMajor;

	if (before == 1) {
		size_t m = activeC;
		size_t k = active; //activeB
		size_t n = after;

		blas::Op op_B = blas::Op::NoTrans;
		blas::gemm(layout, op_h, op_B, m, n, k, alpha, h.coeffs_, m, B.coeffs_, k,
			beta, C.coeffs_, m);

		return;
	}

	size_t m = activeC;
	size_t k = active; // activeB
	size_t n = before;

	size_t pref = before * active;
	size_t prefC = before * activeC;
	blas::Op op_B = blas::Op::Trans;
	T zero = 0.;

	/// Work Memory.
	/// Allocation often dominates computation time, thus,
	/// the memory is static and resize is called to assert
	/// sufficient memory.
	static Tensor<T> mem;
	mem.resize(C.shape_);

	for (size_t aft = 0; aft < after; ++aft) {
		blas::gemm(layout, op_h, op_B, m, n, k, alpha, h.coeffs_, m, &B[aft * pref], n,
			zero, &mem[aft * prefC], m);
		transpose(&C[aft * prefC], &mem[aft * prefC], activeC, before, beta);
	}
}

template void matrixTensor(Tensor<d>& C, const Tensor<d>& h, const Tensor<d>& B,
	size_t k1, d alpha, d beta, blas::Op op_h);
template void matrixTensor(Tensor<cd>& C, const Tensor<cd>& h, const Tensor<cd>& B,
	size_t k1, cd alpha, cd beta, blas::Op op_h);

template<typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& B,
	size_t k, T alpha, T beta, blas::Op op_h) {
	TensorShape shape = B.shape_;
	shape = replaceDimension(shape, k, h.shape_[0]);
	Tensor<T> C(shape);
	matrixTensor(C, h, B, k, alpha, beta, op_h);
	return C;
}

template Tensor<d> matrixTensor(const Tensor<d>& h, const Tensor<d>& B,
	size_t k1, d alpha, d beta, blas::Op op_h);
template Tensor<cd> matrixTensor(const Tensor<cd>& h, const Tensor<cd>& B,
	size_t k1, cd alpha, cd beta, blas::Op op_h);

template<typename T>
void contractionRef(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha, T beta) {
	/// This is a simple implementation of the tensor-hole contraction
	/// which is not optimized for efficiency.
	h *= beta;

	for (size_t aft = 0; aft < bra.shape_.after(k); ++aft) {
		for (size_t actL = 0; actL < bra.shape_[k]; ++actL) {
			for (size_t actR = 0; actR < ket.shape_[k]; ++actR) {
				for (size_t bef = 0; bef < bra.shape_.before(k); ++bef) {
					h(actL, actR) += alpha * bra(bef, actL, aft, k) * ket(bef, actR, aft, k);
				}
			}
		}
	}
}

template void contractionRef(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	size_t k, d alpha, d beta);
template void contractionRef(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t k, cd alpha, cd beta);

template<typename T>
void contraction(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha, T beta) {

	/// Work Memory.
	/// Allocation often dominates computation time, thus,
	/// the memory is static and resize is called to assert
	/// sufficient memory.
	static Tensor<T> bra_mem;
	bra_mem.resize(bra.shape_);
	static Tensor<T> ket_mem;
	ket_mem.resize(ket.shape_);

	size_t A = ket.shape_.before(k);
	size_t B = ket.shape_[k];
	size_t B2 = bra.shape_[k];
	size_t C = ket.shape_.after(k);

	transposeBC(bra_mem.coeffs_, bra.coeffs_, A, B, C);
	transposeBC(ket_mem.coeffs_, ket.coeffs_, A, B2, C);

	size_t AC = A * C;
	blas::Layout layout = blas::Layout::ColMajor;
	blas::Op ct = blas::Op::ConjTrans;
	blas::Op notrans = blas::Op::NoTrans;
	blas::gemm(layout, ct, notrans, B, B2, AC, alpha, bra_mem.coeffs_, AC, ket_mem.coeffs_, AC, beta,
		h.coeffs_, B);
}

template void contraction(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	size_t k, d alpha, d beta);
template void contraction(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t k, cd alpha, cd beta);

template<typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha) {
	Tensor<T> h({bra.shape_[k], ket.shape_[k]});
	contraction(h, bra, ket, k, alpha);
	return h;
}

template Tensor<cd> contraction(const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t k, cd alpha);
template Tensor<d> contraction(const Tensor<d>& bra, const Tensor<d>& ket,
	size_t k, d alpha);

template<typename T>
void qr(Tensor<T>& Q, const Tensor<T>& A) {
	Q = A;
	size_t m = A.shape_.lastBefore();
	size_t n = A.shape_.lastDimension();
	size_t lda = m;
	/// A(m, n);

	/// work memory
	static Tensor<T> tau;
	size_t len = (m > n) ? m : n;
	tau.resize({len});

	CHECK(lapack::geqrf(m, n, Q.coeffs_, lda, tau.coeffs_));

	/// generate Q matrix from reflectors
	CHECK(lapack::ungqr(m, n, n, Q.coeffs_, m, tau.coeffs_));
}

template void qr(Tensor<cd>& Q, const Tensor<cd>& A);
template void qr(Tensor<d>& Q, const Tensor<d>& A);

template<typename T>
Tensor<T> qr(const Tensor<T>& A) {
	Tensor<T> Q;
	qr(Q, A);
	return Q;
}

template Tensor<cd> qr(const Tensor<cd>& A);
template Tensor<d> qr(const Tensor<d>& A);

template<typename T>
Tensor<T> qr(Tensor<T> A, size_t k) {
	Tensor<T> Q(A.shape_);
	transpose(Q, A, k);
	qr(A, Q);
	transpose(Q, A, k, true);
	return Q;
}

template Tensor<d> qr(Tensor<d> A, size_t k);
template Tensor<cd> qr(Tensor<cd> A, size_t k);

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

template void gramSchmidt(Tensor<cd>& A);
template void gramSchmidt(Tensor<d>& A);

template<typename T>
void gramSchmidt(Tensor<T>& A, size_t k) {
	Tensor<T> AT = transpose(A, k);
	gramSchmidt(AT);
	transpose(A, AT, k, true);
}

template void gramSchmidt(Tensor<cd>& A, size_t k);
template void gramSchmidt(Tensor<d>& A, size_t k);

/// SVD
template<typename T>
void svd(Tensor<T>& AthenU, Tensor<T>& VT, Tensord& sigma) {
	const TensorShape& shape = AthenU.shape_;
	size_t m = shape.lastBefore();
	size_t n = shape.lastDimension();
	size_t lda = m;
	size_t ldvt = (m > n) ? n : m; // min(m, n);
	lapack::Job jobu = lapack::Job::OverwriteVec;
	lapack::Job jobvt = lapack::Job::SomeVec;

	lapack::gesvd(jobu, jobvt, m, n, AthenU.coeffs_, lda, sigma.coeffs_, nullptr, m, VT.coeffs_, ldvt);
}

template void svd(Tensor<cd>& AthenU, Tensor<cd>& VT, Tensord& sigma);
template void svd(Tensor<d>& AthenU, Tensor<d>& VT, Tensord& sigma);

template<typename T>
SVD<T> svd(Tensor<T> A) {
	const TensorShape& shape = A.shape_;
	size_t m = shape.lastBefore();
	size_t n = shape.lastDimension();
	n = (m > n) ? n : m;
	Tensor<T> VT({n, n});
	Tensord sigma({n});
	svd(A, VT, sigma);
	return SVD<T>({A, VT, sigma});
}

template SVD<cd> svd(Tensor<cd> A);
template SVD<d> svd(Tensor<d> A);

template<typename T>
SVD<T> svd(const Tensor<T>& A, size_t k) {
	auto AT = transpose(A, k);
	transpose(AT, A, k);
	SVD<T> x = svd(AT);
	get<0>(x) = transpose(get<0>(x), k, true);
	return x;
}

template SVD<cd> svd(const Tensor<cd>& A, size_t k);
template SVD<d> svd(const Tensor<d>& A, size_t k);

/// Eigenvector/values
template<typename T>
void diagonalize(SpectralDecomposition<T>& x) {
	Tensor<T>& A = x.first;
	Tensor<d>& ev = x.second;
	size_t m = A.shape_.lastBefore();
	size_t n = A.shape_.lastDimension();
	if (m != n) {
		cerr << "Error: cannot perform heev on non-diagonal matrix.\n";
		exit(1);
	}

	auto jobz = lapack::Job::Vec;
	auto uplo = lapack::Uplo::Upper;

	CHECK(lapack::heev(jobz, uplo, n, A.coeffs_, n, ev.coeffs_));
}

template void diagonalize(SpectralDecompositioncd& x);
template void diagonalize(SpectralDecompositiond& x);

