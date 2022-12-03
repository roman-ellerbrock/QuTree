//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "benchmark_tensor_holeContraction.h"
#include "benchmark_tensor_matrixTensor.h"
#include "benchmark_tree.h"
//#include "optimize_matrixtensor.h"

namespace benchmark {
	enum benchmark_type {matrixTensor, contraction, tree, peak_gemm};
	enum screening_type {dim, order};

	void run(benchmark_type benchmark, screening_type screening, size_t mode, bool gemm,
	size_t mink, size_t maxk, size_t stepk,
	 mt19937& gen, size_t nsample, ostream& os) {

		if (benchmark == matrixTensor) {
			if (screening == dim) {
				screen_dim_matrixTensor(gen, os, nsample, mode, gemm, mink, maxk, stepk);
			} else if (screening == order) {
				screen_order_matrixTensor(gen, os, nsample, mode, gemm, mink, maxk, stepk);
			} else {
				cerr << "Unknown screening type.\n";
			}
		}

		if (benchmark == contraction) {
			if (screening == dim) {
				screen_dim_contraction(gen, os, nsample, mode, gemm, mink, maxk, stepk);
			} else if (screening == order) {
				screen_order_contraction(gen, os, nsample, mode, gemm, mink, maxk, stepk);
			} else {
				cerr << "Unknown screening type.\n";
			}
		}

		if (benchmark == peak_gemm) {
			benchmark_peak_gemm(gen, os, nsample);
		}

		if (benchmark == tree) {
			screen_nleaves(gen, os, nsample, mink, maxk, stepk);
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
	mt19937 gen(1989);
	size_t benchmark;
	size_t type;
	size_t mode;
	bool gemm;
	size_t nsample;
	size_t mink;
	size_t maxk;
	size_t stepk;

	cout << "Reading parameters...\n";
	cin >> benchmark;
	cin >> type;
	cin >> mode;
	cin >> gemm;
	cin >> nsample;
	cin >> mink;
	cin >> maxk;
	cin >> stepk;

	cout << "Benchmark parameters:\n";
	cout << "benchmark: " << benchmark << endl;
	cout << "type: " << type << endl;
	cout << "mode: " << mode << endl;
	cout << "gemm: " << gemm << endl;
	cout << "nsample: " << nsample << endl;
	cout << "mink: " << mink << endl;
	cout << "maxk: " << maxk << endl;
	cout << "stepk: " << stepk << endl;

	auto benchmark_t = (benchmark::benchmark_type) benchmark;
	auto type_t = (benchmark::screening_type) type;

	benchmark::run(benchmark_t, type_t, mode, gemm, mink, maxk, stepk, gen, nsample, cout);
}

