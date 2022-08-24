//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Tensor/TensorBLAS2.h"
#include <lapack.hh>

typedef float f;
typedef double d;
typedef complex<float> cf;
typedef complex<double> cd;


template<typename T>
void gemm(Tensor<T>& c, const Tensor<T>& a, const Tensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b) {

	/// m: Number of rows of the matrix C and \(op(A)\). m >= 0.
	size_t m = nrows(a.shape_, op_a);
	assert(nrows(c.shape_) == m);

	/// n: Number of cols of the matrix C and \(op(B)\). n >= 0.
	size_t n = ncols(b.shape_, op_b);
	assert(ncols(c.shape_) == n);

	size_t k = ncols(a.shape_, op_a);
	assert(nrows(b.shape_, op_b) == k);

	size_t lda = (op_a == blas::Op::NoTrans) ? m : k;
	size_t ldb = (op_b == blas::Op::NoTrans) ? k : n;
	size_t ldc = m;

	blas::gemm(blas::Layout::ColMajor, op_a, op_b, m, n, k,
		alpha, a.data(), lda, b.data(), ldb, beta, c.data(), ldc);
}

template void gemm(Tensor<cd>& c, const Tensor<cd>& a, const Tensor<cd>& b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b);
template void gemm(Tensor<d>& c, const Tensor<d>& a, const Tensor<d>& b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b);
template void gemm(Tensor<f>& c, const Tensor<f>& a, const Tensor<f>& b,
	f alpha, f beta, blas::Op op_a, blas::Op op_b);

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
Tensor<T> operator*(const Tensor<T>& a, const Tensor<T>& b) {
	return gemm(a, b);
}

template Tensor<f> operator*(const Tensor<f>& a, const Tensor<f>& b);
template Tensor<d> operator*(const Tensor<d>& a, const Tensor<d>& b);
template Tensor<cf> operator*(const Tensor<cf>& a, const Tensor<cf>& b);
template Tensor<cd> operator*(const Tensor<cd>& a, const Tensor<cd>& b);

template<typename T>
Tensor<T> unitarySimilarityTrafo(Tensor<T> a, const Tensor<T>& u) {
	a = gemm(a, u);
	return gemm(adjoint(u), a);
}

template Tensor<cd> unitarySimilarityTrafo(Tensor<cd> a, const Tensor<cd>& u);
template Tensor<d> unitarySimilarityTrafo(Tensor<d> a, const Tensor<d>& u);

template<typename T>
void matrixTensorMode0(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	T alpha, T beta, blas::Op op_h) {

	size_t active = B.shape_[0];
	size_t activeC = C.shape_[0];
	size_t after = B.shape_.after(0);
	size_t m = activeC;
	size_t k = active; //activeB
	size_t n = after;

	blas::Layout layout = blas::Layout::ColMajor;
	blas::Op op_B = blas::Op::NoTrans;
	blas::gemm(layout, op_h, op_B, m, n, k, alpha, h.data(), m, B.data(), k,
		beta, C.data(), m);
}

template void matrixTensorMode0(Tensor<d>& C, const Tensor<d>& h, const Tensor<d>& B,
	d alpha, d beta, blas::Op op_h);
template void matrixTensorMode0(Tensor<cd>& C, const Tensor<cd>& h, const Tensor<cd>& B,
	cd alpha, cd beta, blas::Op op_h);

template<typename T>
void matrixTensorModeD(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	T alpha, T beta, blas::Op op_h) {

	if (op_h == blas::Op::NoTrans) {
		op_h = blas::Op::Trans;
	} else if (op_h == blas::Op::Trans) {
		op_h = blas::Op::NoTrans;
	} else if (op_h == blas::Op::ConjTrans) {
		Matrix<T> h2 = adjoint(h);
		matrixTensorModeD(C, h2, B, alpha, beta, blas::Op::NoTrans);
		return;
	}

	size_t d = B.shape_.lastIdx();
	size_t m = C.shape_.before(d);
	size_t n = C.shape_[d];
	size_t k = B.shape_[d];

	blas::Layout layout = blas::Layout::ColMajor;
	blas::Op op_B = blas::Op::NoTrans;
	blas::gemm(layout, op_B, op_h, 
		m, n, k, 
		alpha, 
		B.data(), m, 
		h.data(), n,
		beta, 
		C.data(), m);
}

template void matrixTensorModeD(Tensor<d>& C, const Tensor<d>& h, const Tensor<d>& B,
	d alpha, d beta, blas::Op op_h);
template void matrixTensorModeD(Tensor<cd>& C, const Tensor<cd>& h, const Tensor<cd>& B,
	cd alpha, cd beta, blas::Op op_h);

template<typename T>
void matrixTensorModeX(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	size_t k1, T alpha, T beta, blas::Op op_h) {

	/// C = alpha * h *(k) B + beta * C
	T beta2 = 0.;
	if (beta == (T) 1.) {
		beta2 = (T) 1.;
	} else if (beta != (T) 0.) {
		cerr << "Cannot call matrixTensorX with beta != {0, 1}\n";
		exit(1);
	}

	blas::Op op_B = blas::Op::NoTrans;
	if (op_h == blas::Op::NoTrans) {
		op_h = blas::Op::Trans;
	} else if (op_h == blas::Op::Trans) {
		op_h = blas::Op::NoTrans;
	} else if (op_h == blas::Op::ConjTrans) {
		Matrix<T> h2 = adjoint(h);
		matrixTensorModeX(C, h2, B, k1, alpha, beta, blas::Op::NoTrans);
		return;
	}

	size_t before = B.shape_.before(k1);
	size_t activeB = B.shape_[k1];
	size_t activeC = C.shape_[k1];
	size_t after = B.shape_.after(k1);

	size_t preB = before * activeB;
	size_t preC = before * activeC;

	size_t m = before;
	size_t n = activeC;
	size_t k = activeB;

	blas::Layout layout = blas::Layout::ColMajor;
	#pragma omp parallel for
	for (size_t aft = 0; aft < after; ++aft) {
		blas::gemm(layout, op_B, op_h, m, n, k, 
			alpha, 
			&B[aft * preB], m, 
			h.data(), n,
			beta2,
			&C[aft * preC], m);
	}
}

template void matrixTensorModeX(Tensor<d>& C, const Tensor<d>& h, const Tensor<d>& B,
	size_t k1, d alpha, d beta, blas::Op op_h);
template void matrixTensorModeX(Tensor<cd>& C, const Tensor<cd>& h, const Tensor<cd>& B,
	size_t k1, cd alpha, cd beta, blas::Op op_h);

template<typename T>
void matrixTensor(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	size_t k1, T alpha, T beta, blas::Op op_h) {
	if (k1 == 0) {
		matrixTensorMode0(C, h, B, alpha, beta, op_h);
	} else if (k1 == B.shape_.lastIdx()) {
		matrixTensorModeD(C, h, B, alpha, beta, op_h);
	} else {
		if (beta == (T) 1. || beta == (T) 0.) {
			matrixTensorModeX(C, h, B, k1, alpha, beta, op_h);
		} else {
			cout << "Warning!\n";
			Tensor<T> mem(C.shape_);
			matrixTensorModeX(mem, h, B, k1, alpha, (T) 0., op_h);
			C *= beta;
			C += mem;
		}
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
	shape.setDimension(h.shape_[0], k);
	Tensor<T> C(shape);
	matrixTensor(C, h, B, k, alpha, beta, op_h);
	return C;
}

template Tensor<d> matrixTensor(const Tensor<d>& h, const Tensor<d>& B,
	size_t k1, d alpha, d beta, blas::Op op_h);
template Tensor<cd> matrixTensor(const Tensor<cd>& h, const Tensor<cd>& B,
	size_t k1, cd alpha, cd beta, blas::Op op_h);


template<typename T>
void contractionMode0(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	T alpha, T beta) {

	// conj(bra) * ket
	size_t BL = bra.shape_[0];
	size_t BR = ket.shape_[0];
	size_t C = ket.shape_.after(0);

	blas::Layout layout = blas::Layout::RowMajor;
	blas::Op ct = blas::Op::ConjTrans;
	blas::Op no = blas::Op::NoTrans;
	Tensor<T> h_add({BR, BL});
	blas::gemm(layout, ct, no, 
	BL, BR, C, 
	alpha, 
	bra.data(), BL,
	ket.data(), BR, 
	beta,
	h_add.data(), BR);
	if (beta == (T) 0.) {
		transpose(h, h_add);
	} else if (beta != (T) 1.) {
		h *= beta;
		h += transpose(h_add);
	} else {
		h += transpose(h_add);
	}
}

template void contractionMode0(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	d alpha, d beta);
template void contractionMode0(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	cd alpha, cd beta);

template<typename T>
void contractionModeD(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	T alpha, T beta) {

	size_t d = ket.shape_.lastIdx();
	size_t BR = ket.shape_[d];
	size_t BL = bra.shape_[d];
	size_t A = ket.shape_.before(d); // C = 1
	blas::Layout layout = blas::Layout::ColMajor;
	blas::Op ct = blas::Op::ConjTrans;
	blas::Op notrans = blas::Op::NoTrans;
	blas::gemm(layout, ct, notrans, 
	BL, BR, A, 
	alpha, 
	bra.data(), A, 
	ket.data(), A, 
	beta,
	h.data(), BL);
}

template void contractionModeD(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	d alpha, d beta);
template void contractionModeD(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	cd alpha, cd beta);


template<typename T>
void contractionModeX(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha, T beta) {
	/**
	 * @brief Tensor hole contraction for an arbitrary index k
	 * 
	 */

	size_t A = ket.shape_.before(k);
	size_t BL = bra.shape_[k];
	size_t BR = ket.shape_[k];
	size_t C = ket.shape_.after(k);

	h *= beta;
	for (size_t c = 0; c < C; ++c) {
		blas::Layout layout = blas::Layout::ColMajor;
		blas::Op ct = blas::Op::ConjTrans;
		blas::Op notrans = blas::Op::NoTrans;
		blas::gemm(layout, ct, notrans, 
		BL, BR, A, 
		alpha, 
		&bra.data()[A * BL * c], A, 
		&ket.data()[A * BR * c], A, 
		1.,
		h.data(), BL);
	}
}

template void contractionModeX(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	size_t k, d alpha, d beta);
template void contractionModeX(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t k, cd alpha, cd beta);


template<typename T>
void contraction(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha, T beta) { 
	/// h = beta * h + alpha * contraction(bra, ket, k)
	/// C = beta * C + alpha * A * C
	if (k == 0) {
		contractionMode0(h, bra, ket, alpha, beta);
	} else if (k == bra.shape_.lastIdx()) {
		contractionModeD(h, bra, ket, alpha, beta);
	} else {
		contractionModeX(h, bra, ket, k, alpha, beta);
	}
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
void contraction(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	const vector<size_t>& holes, T alpha) {
	if (holes.empty()) {
		h = contraction(bra, ket, 0, alpha);
		T val = 0.;
		val = trace(h);
		h = Tensor<T>({1});
		h[0] = (T) val;
	} else if (holes.size() == 1) {
		contraction(h, bra, ket, holes.front(), alpha);
	} else {
		cerr << "contraction only implemented for single or no holes.\n";
		exit(1);
	}
}

template void contraction(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	const vector<size_t>& holes, cd alpha);
template void contraction(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	const vector<size_t>& holes, d alpha);

template<typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket,
	const vector<size_t>& holes, T alpha) {
	if (holes.empty()){
		vector<size_t> lastHole = {bra.shape_.lastIdx()};
		return contraction(bra, ket, lastHole, alpha);
	} else if (holes.size() == 1){
		return contraction(bra, ket, holes.front(), alpha);
	} else {
		cerr << "contraction only implemented for single or no holes.\n";
		exit(1);
	}
}

template Tensor<cd> contraction(const Tensor<cd>& bra, const Tensor<cd>& ket,
	const vector<size_t>& holes, cd alpha);
template Tensor<d> contraction(const Tensor<d>& bra, const Tensor<d>& ket,
	const vector<size_t>& holes, d alpha);

template<typename T>
T dotProduct(const Tensor<T>& bra, const Tensor<T>& ket) {
	return trace(contraction(bra, ket, bra.shape_.lastIdx()));
}

template cd dotProduct(const Tensor<cd>& bra, const Tensor<cd>& ket);
template d dotProduct(const Tensor<d>& bra, const Tensor<d>& ket);
