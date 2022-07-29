#include <cublas.h>
#include <stdafx.h>

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

template <typename TL, typename TR>
__global__ void cudaCast_kernel(TL* dest, const TR* src, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dest[idx] = (TL) src[idx];
    }
}

template <typename TL, typename TR>
void cudaCast(TL* dest, const TR* src, size_t n) {
    constexpr size_t nthreads = 1024;
    size_t nblocks = (n + nthreads - 1) / nthreads;
    cudaCast_kernel<<<nblocks, nthreads>>>(dest, src, n);
}

template void cudaCast(float* dest, const double* src, size_t n);

template <typename T>
__global__ void cudaHadamardProduct_kernel(T* C, const T* A, const T* B, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        C[idx] += A[idx] * B[idx];
    }
}

template <typename T>
void cudaHadamardProduct(T* C, const T* A, const T* B, size_t n) {
    constexpr size_t nthreads = 1024;
    size_t nblocks = (n + nthreads - 1) / nthreads;
    cudaHadamardProduct_kernel<<<nblocks, nthreads>>>(C, A, B, n);
}

template void cudaHadamardProduct(d* C, const d* A, const d* B, size_t n);

template <typename T>
__global__ void cudaDiagMatrixMatrix_kernel(T* C, const T* dA, const T* B, T factor, size_t nrow, size_t ncol) {
//    size_t row = blockIdx.x * blockDim.x + threadIdx.x;
//    size_t col = blockIdx.y * blockDim.y + threadIdx.y;
    size_t id_global = blockIdx.x * blockDim.x + threadIdx.x;
    size_t row = id_global % nrow;
    size_t col = id_global / nrow;
    if (row < nrow && col < ncol) {
        size_t idx = row + nrow * col;
        C[idx] += factor * dA[row] * B[idx];
    }
}

template <typename T>
void cudaDiagMatrixMatrix(T* C, const T* dA, const T* B, T factor, size_t nrow, size_t ncol) {
    constexpr size_t nthreads = 128;
//    size_t nblock_row = (nrow + nthreads - 1) / nthreads;
//    size_t nblock_col = (ncol + nthreads - 1) / nthreads;
//    dim3 nblocks3(nblock_row, nblock_col);
//    dim3 nthreads3(nthreads, nthreads);
//    cout << nrow << " " << ncol << endl;
//    cout << nblocks3.x << " " << nblocks3.y << " " << nblocks3.z << endl;
//    cout << nthreads3.x << " " << nthreads3.y << " " << nthreads3.z << endl;
//    cudaDiagMatrixMatrix_kernel<<<nblocks3, nthreads3>>>(C, dA, B, factor, nrow, ncol);
    size_t nblocksX = (nrow * ncol + nthreads - 1) / nthreads;
    dim3 nblocks3(nblocksX, 1);
    dim3 nthreads3(nthreads, 1);
    cudaDiagMatrixMatrix_kernel<<<nblocks3, nthreads3>>>(C, dA, B, factor, nrow, ncol);
}

template void cudaDiagMatrixMatrix(d* C, const d* dA, const d* B, d factor, size_t nrow, size_t ncol);

/*__global__ void cuda_diagmm(double* C, const double* da, const double* b, size_t nrow, size_t ncol) {
    row=blockIdx.x*blockDim.x+threadIdx.x;
    col=blockIdx.y*blockDim.y+threadIdx.y;
    if (row < nrow && col < ncol) {

    }
}*/

