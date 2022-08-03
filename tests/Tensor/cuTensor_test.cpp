#include <gtest/gtest.h>
#include <iostream>
#include "Util/QMConstants.h"
#include "Tensor/cuTensor.h"
#include "Tensor/TensorBLAS2.h"
#include <chrono>

using namespace std;

double eps = 1e-7;

TEST (cuTensor, Constructor) {
	TensorShape tdim({3, 4, 5, 2});
	cuTensord A(tdim);
		EXPECT_EQ(3 * 4 * 5 * 2, A.shape_.totalDimension());
}

TEST (cuTensor, transfer) {
	using namespace polymorphic;

	Tensord hostA = randomd({ 1000, 1000 });
	cuTensord devA = transfer<double, cuMemory, hostMemory>(hostA);
	Tensord hostA2 = transfer<double, hostMemory, cuMemory>(devA);
	EXPECT_NEAR(0., residual(hostA, hostA2), eps);
}

TEST (cuTensor, zero) {
	Tensord hostA = aranged({ 100, 100 });
	cuTensord devA = transferToGPUd(hostA);
	devA.zero();
	hostA.zero();
	Tensord hostA2 = transferFromGPUd(devA);
	EXPECT_NEAR(0., residual(hostA, hostA2), eps);
}

/*
size_t n_sample = 3;
size_t mini = 1000;
size_t maxi = 2000;
size_t inc_plus = 100;
TEST (pureC, double_dgemmChrono) {
	using namespace chrono;
	time_point<steady_clock> start, end;
 
	cout << "double*:\n";

	double alpha = 1.;
	double beta = 0.;

	for (size_t n = mini; n <= maxi; n += inc_plus) {

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
		double ms = duration_cast<nanoseconds>(end - start).count() / ((double) n_sample * 1000000);
		cout << n << " " << ms << " \n";

		free(A);
		free(B);
		free(C);
	}

}

TEST (cuTensor, dgemmChrono) {
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

	for (size_t n = mini; n <= maxi; n += inc_plus) {
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
		double ms = duration_cast<nanoseconds>(end - start).count() / ((double) n_sample * 1000000);;
		cout << dim << " " << ms << " \n";
	}

}

TEST (Tensor, dgemmChrono) {
	size_t m = 1;
	using namespace chrono;
	time_point<steady_clock> start, end;
 
	cout << "MKL:\n";

	double alpha = 1.;
	double beta = 0.;

	for (size_t n = mini; n <= maxi; n += inc_plus) {
		size_t dim = n * m;
//		Tensor<float> A = random<float>({dim, dim});
//		Tensor<float> B = random<float>({dim, dim});
		Tensor<double> A({dim, dim});
		Tensor<double> B({dim, dim});
		Tensor<double> C({dim, dim});

		start = steady_clock::now();
		for (size_t s = 0; s < n_sample; ++s) {
			gemm(C, A, B, alpha, beta, blas::Op::NoTrans, blas::Op::NoTrans);
		}
		end = steady_clock::now();
		double ms = duration_cast<nanoseconds>(end - start).count() /((double) n_sample * 1000000);
		cout << dim << " " << ms << " \n";
	}
}
*/
