#include "Core/SimultaneousDiagonalization.h"

void SimultaneousDiagonalization::Initialization(vector<FactorMatrixcd>& A, double eps_)
{
	// Number of matrices
	nmat = A.size();
	assert(nmat > 0);

	// Dimension of matrices
	FactorMatrixcd& B = A[0];
	dim = B.Dim1();
	
	// Dimension check
	for (size_t k = 0; k < nmat; k++)
	{
		FactorMatrixcd& B = A[k];
		assert(B.Dim1() == dim);
		assert(B.Dim2() == dim);
	}

	// Set convergence parameter
	eps = eps_;
}

void SimultaneousDiagonalization::Calculate(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo)
{
	bool converged = false;
	int iter = 0;
	int maxiter = 100;

	// Initialize the Diagonalization by rotating to the diagonal
	// representation of one of the matrices in A. This
	// avoids stationary points during the optimization process.
	trafo.Zero();
	for (size_t i = 0; i < dim; i++)
		trafo(i, i) = 1.;

	// This rotates to the eigenspace of the first matrix
	 InitialTransformation(A, trafo);

	// Iterate Jacobirotations until a converged result is reached
	// Measure off-diagonal norm
	double delta = MeasureDiagonality(A);
	double delta_off = MeasureOffDiagonals(A);
//	cout << "Start : " << delta << "\t" << delta_off << endl;
	while (!converged && iter < maxiter)
	{
		// Rotation circle over all elements
		JacobiRotations(A, trafo);

		// Measure off-diagonal norm
		delta = MeasureDiagonality(A);
		delta_off = MeasureOffDiagonals(A);

		// Check convergence
		if (delta < eps) { converged = true; }

//		cout << iter << " : " << delta << "\t" << delta_off << endl;

		iter++;
	}
	if (A.size() > 1)
	{
//		cout << iter << " : " << delta << "\t" << delta_off << endl;
	}
}

void SimultaneousDiagonalization::JacobiRotations(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo)
{
	// Angles for Givens-Rotation
	complex<double> c, s = 0;

	// Swipe over the matrix-dimensions and perform jacobi-rotations
	for (size_t i = 0; i < dim - 1; i++)
	{
		for (size_t j = i + 1; j < dim; j++)
		{
			// Calculate Angles c and s for the elements i and j
			CalculateAngles(c, s, i, j, A);

			assert(abs(1. - abs(c)*abs(c) - abs(s)*abs(s)) < 1E-10);

			// Perform the Givens-Rotation with angles c and s
			RotateMatrices(A, c, s, i, j);

			// Rotate the Transformation-Matrix
			GivensTrafoRotation(trafo, c, s, i, j);
		}
	}
}

double SimultaneousDiagonalization::MeasureOffDiagonals(const vector<FactorMatrixcd>& A)
{
	// Measure the norm of off-diagonal elements
	double eps = 0;

	for (size_t k = 0; k < A.size(); k++)
	{
		FactorMatrixcd B(A[k]);

		for (int n = 0; n < B.Dim(); n++)
			B(n, n) = 0;

		B = B*B;
		eps += real(B.Trace());
	}

	return sqrt(eps);
}

double SimultaneousDiagonalization::MeasureDiagonality(vector<FactorMatrixcd>& A)
{
	// Measure the norm of off-diagonal elements
	double eps = 0;

	for (size_t k = 0; k < nmat; k++)
	{
		FactorMatrixcd&B = A[k];

		for (size_t n = 0; n < dim; n++)
			for (size_t m = 0; m < dim; m++)
				if (m != n)
				{
					eps += pow(abs(B(m, n)), 2);
				}
	}

	return sqrt(eps);
}

void SimultaneousDiagonalization::InitialTransformation(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo)
{
	FactorMatrixcd& B = A[0];
	Vectord ev(B.Dim1());
	B.cDiag(trafo, ev);

	for (size_t k = 0; k < A.size(); k++)
	{
		FactorMatrixcd& C = A[k];
		int mode = C.Mode();
		Matrixcd Cmat= UnitarySimilarityTrafo(C, trafo);
		C = FactorMatrixcd(Cmat, mode);
	}
}
