#include "cuTensor.h"
#include "TensorBLAS1.h"
#include <lapack.hh>

typedef float f;
typedef double d;
typedef complex<double> cd;

//using namespace polymorphic;
//template class Tensor<float, cuMemory>;
//template class Tensor<double, cuMemory>;
//template class Tensor<complex<double>, cuMemory>;

using namespace polymorphic;

template<typename T>
void gemm(cuTensor<T>& c, const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b, blas::Queue& queue) {

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
		alpha, a.data(), lda, b.data(), ldb, beta, c.data(), ldc,
		queue);
}

template void gemm(cuTensor<cd>& c, const cuTensor<cd>& a, const cuTensor<cd>& b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template void gemm(cuTensor<d>& c, const cuTensor<d>& a, const cuTensor<d>& b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template void gemm(cuTensor<f>& c, const cuTensor<f>& a, const cuTensor<f>& b,
	f alpha, f beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);


template<typename T, template <typename> class Mem, template <typename> class oMem>
Tensor<T, Mem> transfer(const Tensor<T, oMem>& src) {
	Tensor<T, Mem> dest(src.shape_);
	dest.mem() = src.mem();
	return dest;
}

template Tensor<d, hostMemory> transfer(const Tensor<d, cuMemory>& src);
template Tensor<cd, hostMemory> transfer(const Tensor<cd, cuMemory>& src);
template Tensor<d, cuMemory> transfer(const Tensor<d, hostMemory>& src);
template Tensor<cd, cuMemory> transfer(const Tensor<cd, hostMemory>& src);
