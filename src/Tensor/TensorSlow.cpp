//
// Created by Roman Ellerbrock on 11/17/21.
//
#include "Tensor/TensorSlow.h"
#include "Tensor/TensorBLAS1.h"

typedef complex<double> cd;
typedef double d;

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
	assert(c.shape_[0] == a.shape_[0]);
	assert(c.shape_[1] == b.shape_[1]);
	assert(a.shape_[1] == b.shape_[0]);
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

template<typename T>
Tensor<T> gemmRef(Tensor<T> a, Tensor<T> b,
	T alpha, blas::Op op_a, blas::Op op_b) {

	size_t m = nrows(a.shape_, op_a);
	size_t n = ncols(b.shape_, op_b);
	Tensor<T> c({m, n});
	T beta = 0.;
	gemmRef(c, a, b, alpha, beta, op_a, op_b);
	return c;
}

template Tensor<cd> gemmRef(Tensor<cd> a, Tensor<cd> b,
	cd alpha, blas::Op op_a, blas::Op op_b);
template Tensor<d> gemmRef(Tensor<d> a, Tensor<d> b,
	d alpha, blas::Op op_a, blas::Op op_b);

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

	/// C[bef, actC, aft] += h[actC, act] * B[bef, act, aft]
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
void contractionRef(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha, T beta) {
	/// This is a simple implementation of the tensor-hole contraction
	/// which is not optimized for efficiency.
	if (bra.shape_[k] != h.shape_[0]) {
		cerr << "Wrong dimension in contractionRef (l).\n";
	}
	if (ket.shape_[k] != h.shape_[1]) {
		cerr << "Wrong dimension in contractionRef (r).\n";
	}
	h *= beta;

	for (size_t aft = 0; aft < bra.shape_.after(k); ++aft) {
		for (size_t actL = 0; actL < bra.shape_[k]; ++actL) {
			for (size_t actR = 0; actR < ket.shape_[k]; ++actR) {
				for (size_t bef = 0; bef < bra.shape_.before(k); ++bef) {
					h(actL, actR) += alpha * conj(bra(bef, actL, aft, k)) * ket(bef, actR, aft, k);
				}
			}
		}
	}
}

template void contractionRef(Tensor<d>& h, const Tensor<d>& bra, const Tensor<d>& ket,
	size_t k, d alpha, d beta);
template void contractionRef(Tensor<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t k, cd alpha, cd beta);

template <typename T>
Tensor<T> toTensorRef(const SVD<T>& svd) {
	const Tensor<T>& U = svd.U();
	Tensor<T> VT = svd.VT();
	const Tensord& sigma = svd.sigma();

	for (size_t i = 0; i < VT.shape_[0]; ++i) {
		for (size_t j = 0; j < VT.shape_[1]; ++j) {
			VT(i, j) *= sigma(i);
		}
	}
	return gemmRef(U, VT);
}

template Tensor<cd> toTensorRef(const SVD<cd>& svd);
template Tensor<d> toTensorRef(const SVD<d>& svd);

template<typename T>
Tensor<T> toTensorRef(const SpectralDecomposition<T>& X) {
	const Tensor<T>& mat = X.U();
	const Tensord& vec = X.ev();
	auto mat2(mat);
	for (size_t i = 0; i < mat2.shape_[0]; ++i) {
		for (size_t k = 0; k < mat2.shape_[1]; ++k) {
			mat2(i, k) *= vec(k);
		}
	}
	return gemmRef(mat, adjoint(mat2));
}

template Tensor<cd> toTensorRef(const SpectralDecomposition<cd>& X);
template Tensor<d> toTensorRef(const SpectralDecomposition<d>& X);
