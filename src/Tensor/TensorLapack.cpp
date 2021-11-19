//
// Created by Roman Ellerbrock on 11/18/21.
//
#include "Tensor/TensorLapack.h"
#include "Check.h"

typedef complex<double> cd;

typedef double d;

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
Tensor<T> toTensor(const SVD<T>& x) {
	const auto& U = x.U();
	auto VT = x.VT();
	const auto& sigma = x.sigma();
	assert(U.shape_.order() == 2);
	assert(U.shape_[1] == sigma.shape_[0]);
	assert(VT.shape_[1] == sigma.shape_[0]);
	for (size_t i = 0; i < VT.shape_[0]; ++i) {
		for (size_t j = 0; j < VT.shape_[1]; ++j) {
			VT(i, j) *= sigma(i);
		}
	}
	return gemm(U, VT);
}

template Tensor<cd> toTensor(const SVD<cd>& x);
template Tensor<d> toTensor(const SVD<d>& x);

template<typename T>
void svd(SVD<T>& x, Tensor<T> A) {
	Tensor<T>& U = x.U();
	Tensor<T>& VT = x.VT();
	Tensord& sigma = x.sigma();

	const TensorShape& shape = A.shape_;
	size_t m = shape.lastBefore();
	size_t n = shape.lastDimension();
	size_t lda = m;
	size_t ldvt = (m > n) ? n : m; // min(m, n);
	size_t ldu = U.shape_.lastBefore();

	lapack::Job jobu = lapack::Job::SomeVec;
	lapack::Job jobvt = lapack::Job::SomeVec;

	lapack::gesvd(jobu, jobvt, m, n, A.coeffs_, lda, sigma.coeffs_, U.coeffs_, ldu, VT.coeffs_, ldvt);
}

template void svd(SVDcd&, Tensorcd);
template void svd(SVDd&, Tensord);

template<typename T>
SVD<T> svd(Tensor<T> A) {
	SVD<T> x(A.shape_);
	svd(x, A);
	return x;
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
Tensor<T> toTensor(const SpectralDecomposition<T>& X) {
	const Tensor<T>& mat = X.U();
	const Tensor<d>& vec = X.ev();

	auto mat2(mat);
	for (size_t i = 0; i < mat2.shape_[0]; ++i) {
		for (size_t k = 0; k < mat2.shape_[1]; ++k) {
			mat2(i, k) *= vec(k);
		}
	}
	return gemm(mat, adjoint(mat2));
}

template Tensor<cd> toTensor(const SpectralDecomposition<cd>& X);
template Tensor<d> toTensor(const SpectralDecomposition<d>& X);

template<typename T>
void heev(SpectralDecomposition<T>& x) {
	Tensor<T>& A = x.U();
	Tensor<d>& ev = x.ev();
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

template void heev(SpectralDecompositioncd& x);
template void heev(SpectralDecompositiond& x);

