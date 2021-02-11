#include "Util/SimultaneousDiagonalization.h"

void SimultaneousDiagonalization::Initialization(vector<Matrixcd>& A,
	double eps_) {
	// Number of matrices
	nmat_ = A.size();
	assert(nmat_ > 0);

	// Dimension of matrices
	Matrixcd& B = A[0];
    dim_ = B.dim1();

	// Dimension check
	for (size_t k = 0; k < nmat_; k++) {
		Matrixcd& B = A[k];
		assert(B.dim1() == dim_);
		assert(B.dim2() == dim_);
	}

	// Set convergence parameter
	eps_ = eps_;
}

void SimultaneousDiagonalization::Calculate(vector<Matrixcd>& A,
	Matrixcd& trafo) {
	bool converged = false;
	int iter = 0;
	int maxiter = 100;

	// Initialize the Diagonalization by rotating to the diagonal
	// representation of one of the matrices in A. This
	// avoids stationary points during the optimization process.
	trafo.zero();
	for (size_t i = 0; i < dim_; i++)
		trafo(i, i) = 1.;

	// This rotates to the eigenspace of the first matrix
	InitialTransformation(A, trafo);

	// Iterate Jacobirotations until a converged result is reached
	// Measure off-diagonal norm
	double delta = MeasureDiagonality(A);
	double delta_off = MeasureOffDiagonals(A);
//	cout << "Start : " << delta << "\t" << delta_off << endl;
	while (!converged && iter < maxiter) {
		// Rotation circle over all elements
		JacobiRotations(A, trafo);

		// Measure off-diagonal norm
		delta = MeasureDiagonality(A);
		delta_off = MeasureOffDiagonals(A);

		// Check convergence
		if (delta < eps_) { converged = true; }

//		cout << iter << " : " << delta << "\t" << delta_off << endl;

		iter++;
	}
	if (A.size() > 1) {
//		cout << iter << " : " << delta << "\t" << delta_off << endl;
	}
}

void SimultaneousDiagonalization::JacobiRotations(vector<Matrixcd>& A,
	Matrixcd& trafo) {
	// Angles for Givens-Rotation
	complex<double> c, s = 0;

	// Swipe over the matrix-dimensions and perform jacobi-rotations
	for (size_t i = 0; i < dim_ - 1; i++) {
		for (size_t j = i + 1; j < dim_; j++) {
			// Calculate Angles c and s for the elements i and j
			CalculateAngles(c, s, i, j, A);

			assert(abs(1. - abs(c) * abs(c) - abs(s) * abs(s)) < 1E-10);

			// Perform the Givens-Rotation with angles c and s
			RotateMatrices(A, c, s, i, j);

			// Rotate the Transformation-Matrix
			GivensTrafoRotation(trafo, c, s, i, j);
		}
	}
}

double SimultaneousDiagonalization::MeasureOffDiagonals(const vector<Matrixcd>& A) {
	// Measure the norm of off-diagonal elements
	double eps = 0;

	for (size_t k = 0; k < A.size(); k++) {
		Matrixcd B(A[k]);

		for (int n = 0; n < B.dim1(); n++)
			B(n, n) = 0;

		B = B * B;
		eps += real(B.trace());
	}

	return sqrt(eps);
}

double SimultaneousDiagonalization::MeasureDiagonality(vector<Matrixcd>& A) {
	// Measure the norm of off-diagonal elements
	double eps = 0;

	for (size_t k = 0; k < nmat_; k++) {
		Matrixcd& B = A[k];

		for (size_t n = 0; n < dim_; n++)
			for (size_t m = 0; m < dim_; m++)
				if (m != n) {
					eps += pow(abs(B(m, n)), 2);
				}
	}

	return sqrt(eps);
}

void SimultaneousDiagonalization::InitialTransformation(vector<Matrixcd>& A,
	Matrixcd& trafo) {
	Matrixcd& B = A[0];
	Vectord ev(B.dim1());
	B.cDiag(trafo, ev);

	for (size_t k = 0; k < A.size(); k++) {
		A[k] = unitarySimilarityTrafo(A[k], trafo);
	}
}
