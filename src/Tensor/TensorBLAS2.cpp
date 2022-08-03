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
		blas::gemm(layout, op_h, op_B, m, n, k, alpha, h.data(), m, B.data(), k,
			beta, C.data(), m);

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
		blas::gemm(layout, op_h, op_B, m, n, k, alpha, h.data(), m, &B[aft * pref], n,
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

	transposeBC(bra_mem.data(), bra.data(), A, B2, C);
	transposeBC(ket_mem.data(), ket.data(), A, B, C);

	size_t AC = A * C;
	blas::Layout layout = blas::Layout::ColMajor;
	blas::Op ct = blas::Op::ConjTrans;
	blas::Op notrans = blas::Op::NoTrans;
	blas::gemm(layout, ct, notrans, 
	B2, B, AC, 
	alpha, 
	bra_mem.data(), AC, 
	ket_mem.data(), AC, 
	beta,
	h.data(), B2);
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
