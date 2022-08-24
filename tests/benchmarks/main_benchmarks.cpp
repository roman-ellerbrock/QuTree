//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "benchmark_tensor.h"
//#include "optimize_matrixtensor.h"

namespace benchmark {
	enum benchmark_type {matrixTensor, contraction, tree, peak_gemm};
	enum screening_type {dim, order};

	void run() {
		mt19937 gen(1989);
		size_t nsample = 4;
		ostream& os = cout;

		auto benchmark = peak_gemm;
		auto screening = order;
		size_t mode = 2;
		bool gemm = true;

		if (benchmark == matrixTensor) {
			if (screening == dim) {
				screen_dim_matrixTensor(gen, os, nsample, mode, gemm);
			} else if (screening == order) {
				screen_order_matrixTensor(gen, os, nsample, mode, gemm);
			} else {
				cerr << "Unknown screening type.\n";
			}
		}

		if (benchmark == contraction) {
			if (screening == dim) {
				screen_dim_contraction(gen, os, nsample, mode, gemm);
			} else if (screening == order) {
				screen_order_contraction(gen, os, nsample, mode, gemm);
			} else {
				cerr << "Unknown screening type.\n";
			}
		}

		if (benchmark == peak_gemm) {
			benchmark_peak_gemm(gen, os, nsample);
		}

//		screen_order(gen, os, nsample);
//		screen_dim(gen, os, nsample);
//		benchmark_peak_gemm(gen, os, nsample);
//		screen_nleaves(gen, os, nsample);
//		cout << "Transpose:\n";
//		screenTranspose(gen, os, nsample);
//		cout << "TransposeAB:\n";
//		screenDimensionTransposeAB(gen, os, nsample);
//		cout << "matrixTensor:\n";
//		screenDimensionMatrixTensor(gen, os, nsample);
//		cout << "TensorContraction:\n";
//		screenDimensionTensorHoleProduct(gen, os, nsample);
//		cout << "DGEEM:\n";
//		screenDGEEM(gen, os, nsample);
	}
}

int main() {
	benchmark::run();
}

