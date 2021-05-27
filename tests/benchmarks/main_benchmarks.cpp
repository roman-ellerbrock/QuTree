//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "benchmark_tensor.h"
#include "optimize_matrixtensor.h"

namespace benchmark {
	void run() {
		mt19937 gen(1989);
		size_t nsample = 20;
		ostream& os = cout;

//		screen_order(gen, os, nsample);
//		screen_dim(gen, os, nsample);
//		screen_nleaves(gen, os, nsample);
//		cout << "Transpose:\n";
//		screenTranspose(gen, os, nsample);
//		cout << "TransposeAB:\n";
//		screenDimensionTransposeAB(gen, os, nsample);
//		cout << "matrixTensor:\n";
//		screenDimensionMatrixTensor(gen, os, nsample);
//		cout << "TensorContraction:\n";
		screenDimensionTensorHoleProduct(gen, os, nsample);
	}
}

int main() {
	benchmark::run();
}

