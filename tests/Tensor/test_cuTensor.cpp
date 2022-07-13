#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Util/QMConstants.h"
#include "Tensor/cuTensor.h"
#include "Tensor/TensorBLAS2.h"
#include <chrono>

using namespace std;

SUITE (cuTensor) {

	double eps = 1e-7;

	TEST (cuTensor_Constructor) {
		TensorShape tdim({3, 4, 5, 2});
		cuTensord A(tdim);
			CHECK_EQUAL(3 * 4 * 5 * 2, A.shape_.totalDimension());
	}

	TEST (cuTensor_dgemm) {
		cuTensord A({100, 100});
		cuTensord B({100, 100});
		cuTensord C({100, 100});

		int device = 0;	
		int batch_size = 1000;
		blas::Queue queue(device, batch_size);
		blas::set_device(device);
		gemm(C, A, B, 1., 0., blas::Op::NoTrans, blas::Op::NoTrans, queue);
	}

	TEST (Tensor_transfer) {
		using namespace polymorphic;

		Tensord hostA = randomd({ 1000, 1000 });
		cuTensord devA = transfer<double, cuMemory, hostMemory>(hostA);
		Tensord hostA2 = transfer<double, hostMemory, cuMemory>(devA);
		CHECK_CLOSE(0., residual(hostA, hostA2), eps);
	}

	TEST (overhead) {
		cuTensord A({100, 100});
		cuTensord B({100, 100});
		cuTensord C({100, 100});

		int device = 0;	
		int batch_size = 1000;
		blas::Queue queue(device, batch_size);
		blas::set_device(device);
		gemm(C, A, B, 1., 0., blas::Op::NoTrans, blas::Op::NoTrans, queue);
		
	}

	size_t n_sample = 10;
	size_t mini = 1000;
	size_t maxi = 2000;
	TEST (double_dgemmChrono) {
		using namespace chrono;
		time_point<steady_clock> start, end;
  
		cout << "double*:\n";

		double alpha = 1.;
		double beta = 0.;

		for (size_t n = mini; n <= maxi; n += 1000) {
			double* A = (double*) malloc(n* n * sizeof(double));
			double* B = (double*) malloc(n* n * sizeof(double));
			double* C = (double*) malloc(n* n * sizeof(double));

			start = steady_clock::now();
			for (size_t s = 0; s < n_sample; ++s) {
			using namespace blas;
				gemm(Layout::ColMajor, Op::NoTrans, Op::NoTrans,
					n, n, n, alpha, A, n, B, n, beta, C, n);
			}
			end = steady_clock::now();
			double ms = duration_cast<nanoseconds>(end - start).count() / 1000000;
			cout << n << " " << ms << " \n";
			free(A);
			free(B);
			free(C);
		}

	}

	TEST (cuTensor_dgemmChrono) {
		size_t m = 1;
		using namespace chrono;
		time_point<steady_clock> start, end;
  
		cout << "CUDA:\n";

		int device = 0;	
		int batch_size = 2048;
		blas::Queue queue(device, batch_size);
		blas::set_device(device);

		double alpha = 1.;
		double beta = 0.;

		for (size_t n = mini; n <= maxi; n += 1000) {
			size_t dim = n * m;
			cuTensor<double> A({dim, dim});
			cuTensor<double> B({dim, dim});
			cuTensor<double> C({dim, dim});

			queue.sync();
			start = steady_clock::now();
			for (size_t s = 0; s < n_sample; ++s) {
				gemm(C, A, B, alpha, beta, blas::Op::NoTrans, blas::Op::NoTrans, queue);
			}
			queue.sync();
			end = steady_clock::now();
			double ms = duration_cast<nanoseconds>(end - start).count() / 1000000;
			cout << dim << " " << ms << " \n";
		}

	}

	TEST (Tensor_dgemmChrono) {
		size_t m = 1;
		using namespace chrono;
		time_point<steady_clock> start, end;
  
		cout << "MKL:\n";

		double alpha = 1.;
		double beta = 0.;

		for (size_t n = mini; n <= maxi; n += 1000) {
			size_t dim = n * m;
//			Tensor<float> A = random<float>({dim, dim});
//			Tensor<float> B = random<float>({dim, dim});
			Tensor<double> A({dim, dim});
			Tensor<double> B({dim, dim});
			Tensor<double> C({dim, dim});

			start = steady_clock::now();
			for (size_t s = 0; s < n_sample; ++s) {
				gemm(C, A, B, alpha, beta, blas::Op::NoTrans, blas::Op::NoTrans);
			}
			end = steady_clock::now();
			double ms = duration_cast<nanoseconds>(end - start).count() / 1000000;
			cout << dim << " " << ms << " \n";
		}
	}


}
