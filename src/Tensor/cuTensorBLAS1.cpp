#include "Tensor/TensorBLAS1.hpp"
#include "Tensor/cuTensorBLAS1.h"
#include "Tensor/cuTensorBLAS1.cuh"

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

//template double nrm2<d, cuTensor, blas::Queue>(const cuTensor<d>& A, size_t inc_a, blas::Queue&);

template void axpy<f, cuTensor<f>, blas::Queue>(const cuTensor<f>& A, cuTensor<f>& B, f alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<d, cuTensor<d>, blas::Queue>(const cuTensor<d>& A, cuTensor<d>& B, d alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cf, cuTensor<cf>, blas::Queue>(const cuTensor<cf>& A, cuTensor<cf>& B, cf alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cd, cuTensor<cd>, blas::Queue>(const cuTensor<cd>& A, cuTensor<cd>& B, cd alpha, size_t inc_a, size_t inc_b, blas::Queue&);

/*template <>
void cast<float, double, polymorphic::cuMemory>(cuTensord& L, const cuTensord& R) {
    size_t n = L.shape_.totalDimension();
    cuCastdf(L.data(), R.data(), n);
}
*/

//template <typename, typename,template <typename> class>
template <>
void cast<float, double, polymorphic::cuMemory>(cuTensorf& L, const cuTensord& R) {
    cudaCast(L.data(), R.data(), L.shape_.totalDimension());
}
