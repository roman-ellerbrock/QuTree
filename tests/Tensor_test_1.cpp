#include "Tensor.h"
#include <iostream>

using namespace std;

int test_TensorDim_1() {
	/// Create a TensorDim, write to file, read in again
	TensorDim tdim({3, 4, 5}, 2);
	ofstream file("tdim.dat");
	tdim.WriteBin(file);
	file.close();
	TensorDim odim("tdim.dat");
	return (odim == tdim);
}

bool test_TensorDim_2() {
	/// Check Getters and Initialization
	bool success = true;

	TensorDim tdim({3, 4, 5}, 2);

	if (tdim.getdimtot() != 3 * 4 * 5 * 2) { success = false; }
	if (tdim.getdimpart() != 3 * 4 * 5) { success = false; }
	if (tdim.getntensor() != 2) { success = false; }

	return success;
}

bool test_Matrix_1() {
	/// Test Matrix I/O
	size_t dim1 = 2;
	size_t dim2 = 3;
	Matrixcd M(dim1, dim2);
	for (size_t i = 0; i < dim1; ++i) {
		for (size_t j = 0; j < dim2; ++j) {
			M(i, j) = j * i;
		}
	}
	M.Write("matrix1.dat");
	Matrixcd N("matrix1.dat");
	return (M == N);
}

bool test_Tensor_1() {
	/// Test Tensor I/O
	bool success = true;

	TensorDim tdim({3, 4, 5}, 2);
	Tensorcd A(tdim);
	for (size_t i = 0; i < tdim.getdimtot(); ++i) {
		A(i) = i;
	}
	ofstream file("tensor1.dat");
	A.Write(file);
	file.close();
	Tensorcd B("tensor1.dat");
	auto C = A - B;
	auto s = C.DotProduct(C);
	auto norm = abs(s.Trace());
	if (norm > 1e-9) {
		success = false;
		cout << "Failed non-binary I/O test.\n";
	}
	return success;
}

Tensorcd test_Tensor_0() {
	TensorDim tdim({3, 4, 5}, 2);
	Tensorcd A(tdim);
	return A;
}

void example_Tensor_0() {
	TensorDim tdim({3, 4, 5});
	Tensorcd A(tdim);
	Tensorcd B(tdim);
	auto B0 = HoleProduct(A, B, 0);

	auto spectral_Decomp_0 = B0.cDiag();
	auto eigenvectors0 = spectral_Decomp_0.second;
	auto U0 = spectral_Decomp_0.first;

	TensorDim tdim2({3, 20});
	B.Reshape(tdim2);

	TensorDim dim({3, 4, 5, 6});
	Tensorcd X(dim);
	Matrixcd s(4, 4);
	auto Y = multAB(s, X, 1);
}

void testi() {
	auto C = test_Tensor_0(); // Call copy asignment operator=()
	auto D(C); // Call copy constructor
	D[0] = 1.;
	auto x = D[1];
}

int test_Tensor() {
	cout << "Testing Matrix:\n";
	cout << "I/O:\n";
	cout << test_Matrix_1() << endl;

	cout << "Testing TensorDim:\n";
	cout << "I/O:\n";
	cout << test_TensorDim_1() << endl;
	cout << "Getters:\n";
	cout << test_TensorDim_2() << endl;

	cout << "Testing Tensor I/O\n";
	cout << test_Tensor_1() << endl;

	return 1;
}


