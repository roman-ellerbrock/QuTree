#include <cublas.h>
#include <stdafx.h>

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

/// ================================================================================ 

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
template void cudaCast(double* dest, const float* src, size_t n);

/// ================================================================================ 

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

/// ================================================================================ 

template <typename T>
__global__ void cudaDiagMatrixMatrix_kernel(T* C, const T* dA, const T* B, T factor, size_t nrow, size_t ncol) {
    size_t row = blockIdx.x * blockDim.x + threadIdx.x;
    size_t col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row < nrow && col < ncol) {
        size_t idx = row + nrow * col;
        C[idx] += factor * dA[row] * B[idx];
    }
}

template <typename T>
void cudaDiagMatrixMatrix(T* C, const T* dA, const T* B, T factor, size_t nrow, size_t ncol) {
    constexpr size_t nthreads = 32;
    size_t nblocksX = (nrow + nthreads - 1) / nthreads;
    size_t nblocksY = (ncol + nthreads - 1) / nthreads;
    dim3 nblocks3(nblocksX, nblocksY);
    dim3 nthreads3(nthreads, nthreads);
    cudaDiagMatrixMatrix_kernel<<<nblocks3, nthreads3>>>(C, dA, B, factor, nrow, ncol);
}

template void cudaDiagMatrixMatrix(d* C, const d* dA, const d* B, d factor, size_t nrow, size_t ncol);

/// ================================================================================ 

template <typename T>
__global__ void cudaMatrixDiagMatrix_kernel(T* C, const T* A, const T* dB, T factor, size_t nrow, size_t ncol) {
    size_t row = blockIdx.x * blockDim.x + threadIdx.x;
    size_t col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row < nrow && col < ncol) {
        size_t idx = row + nrow * col;
        C[idx] += factor * A[idx] * dB[col];
    }
}

template <typename T>
void cudaMatrixDiagMatrix(T* C, const T* A, const T* dB, T factor, size_t nrow, size_t ncol) {
    constexpr size_t nthreads = 32;
    size_t nblocksX = (nrow + nthreads - 1) / nthreads;
    size_t nblocksY = (ncol + nthreads - 1) / nthreads;
    dim3 nblocks3(nblocksX, nblocksY);
    dim3 nthreads3(nthreads, nthreads);
    cudaMatrixDiagMatrix_kernel<<<nblocks3, nthreads3>>>(C, A, dB, factor, nrow, ncol);
}

template void cudaMatrixDiagMatrix(d* C, const d* A, const d* dB, d factor, size_t nrow, size_t ncol);

/// ================================================================================ 

/// C += factor * (A * d + d * A)
template <typename T>
__global__ void cudaDiagmPlusmdiag_kernel(T* C, const T* A, const T* diag, T factor, size_t nrow, size_t ncol) {
    size_t row = blockIdx.x * blockDim.x + threadIdx.x;
    size_t col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row < nrow && col < ncol) {
        size_t idx = row + nrow * col;
        C[idx] += factor * A[idx] * (diag[col] + diag[row]);
    }
}

template <typename T>
void cudaDiagmPlusmdiag(T* C, const T* A, const T* diag, T factor, size_t nrow, size_t ncol) {
    constexpr size_t nthreads = 32;
    size_t nblocksX = (nrow + nthreads - 1) / nthreads;
    size_t nblocksY = (ncol + nthreads - 1) / nthreads;
    dim3 nblocks3(nblocksX, nblocksY);
    dim3 nthreads3(nthreads, nthreads);
    cudaDiagmPlusmdiag_kernel<<<nblocks3, nthreads3>>>(C, A, diag, factor, nrow, ncol);
}

template void cudaDiagmPlusmdiag(d* C, const d* A, const d* diag, d factor, size_t nrow, size_t ncol);


/*__global__ void cuda_diagmm(double* C, const double* da, const double* b, size_t nrow, size_t ncol) {
    row=blockIdx.x*blockDim.x+threadIdx.x;
    col=blockIdx.y*blockDim.y+threadIdx.y;
    if (row < nrow && col < ncol) {

    }
}*/

