#include <cublas.h>
#include <stdafx.h>

template <typename TL, typename TR>
__global__ void cudaCast_kernel(TL* dest, const TR* src, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dest[idx] = (float) src[idx];
    }
}

template <typename TL, typename TR>
void cudaCast(TL* dest, const TR* src, size_t n) {
    constexpr size_t nthreads = 1024;
    size_t nblocks = (n + nthreads - 1) / nthreads;
    cudaCast_kernel<<<nblocks, nthreads>>>(dest, src, n);
}

template void cudaCast(float* dest, const double* src, size_t n);

/*__global__ void cuda_diagmm(double* C, const double* da, const double* b, size_t nrow, size_t ncol) {
    row=blockIdx.x*blockDim.x+threadIdx.x;
    col=blockIdx.y*blockDim.y+threadIdx.y;
    if (row < nrow && col < ncol) {

    }
}*/