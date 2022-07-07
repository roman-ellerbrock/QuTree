#include "cblas.h"
#include "lapacke.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cublas_v2.h>
#include <curand.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void simple_dgemm(float* C, float* A, float* B, int m) {
    float alpha = 1.0f;
    float beta = 0.0f;
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, alpha, A, m, B, m, beta, C, m);
}

/*void simple_cublasDgemm(float* C, float* A, float* B, int m) {
        float alpha = 1.;
    float beta = 0.;
    cublasDgemm(cublasGlobalHandle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m,
        (float *)gpu_alpha1, (float *)&gpu_A[0], 1024, (float *)&gpu_B
         [0], 1024, (float *)gpu_beta1, (float *)&gpu_C[0], 1024);
}*/

float residual(float* A, float* B, int N) {
    float r = 0.;
    for (int i = 0; i < N; ++i) {
        r += pow(abs(A[i] - B[i]), 2);
    }
    return sqrt(r);
}

int cu_main(int m, int runs) {
    int N = m*m;
    int bytes = N * sizeof(float);

    float *dev_A, *dev_B, *dev_C;
    float* host_A = (float*) malloc(bytes);
    float* host_B = (float*) malloc(bytes);
    float* host_C = (float*) malloc(bytes);
    float* host_C2 = (float*) malloc(bytes);
    cudaMalloc(&dev_A, bytes);
    cudaMalloc(&dev_B, bytes);
    cudaMalloc(&dev_C, bytes);

    // random number generator
    curandGenerator_t prng;
    curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);

    // Set the seed
    curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock());

    // Fill matrices with random numbers on device
    curandGenerateUniform(prng, dev_A, N);
    curandGenerateUniform(prng, dev_B, N);

    // cuBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);

    float alpha = 1.0f;
    float beta = 0.0f;
    
    // prepare CUDA event timer
    cudaEvent_t start, stop;
    float time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // call dgemm on Device
    cudaEventRecord(start, 0);
    for (size_t i = 0; i < runs; ++i) {
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &alpha, dev_A, m, dev_B, m, &beta, dev_C, m);
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    gpuErrchk(cudaPeekAtLastError());

    // report time
    cudaEventElapsedTime(&time, start, stop);
    double d_t = time/= (double) runs;
//    printf ("Time on device: %f ms\n", d_t);

    // copy data back and compare
    cudaMemcpy(host_A, dev_A, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_B, dev_B, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_C, dev_C, bytes, cudaMemcpyDeviceToHost);
    
    // prepare CPU timer
    clock_t t;
    t = clock();

    // cal dgemm on Host
    for (size_t i = 0; i < runs; ++i) {
        simple_dgemm(host_C2, host_A, host_B, m);
    }
    t = clock() - t;
    double h_t = ((float)t*1000)/CLOCKS_PER_SEC/(double)runs;
//    printf ("Time on host: %f ms.\n", h_t);

    // verify solution
//    std::cout << "Residual: " << residual(host_C, host_C2, N) << std::endl;

    std::cout << m << " " << h_t << " " << d_t << std::endl;
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_C);
    return 0;
}

 