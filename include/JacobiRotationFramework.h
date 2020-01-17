#pragma once
#include "SingleParticleOperator.h"

namespace JacobiRotationFramework
{
	// Perform a Givens-rotation on an SPO 
	void GivensRotation(SPOcd& B, complex<double> c, complex<double> s,
	                    int i, int j);

	// Perform a Givens-Rotation on a set of SPOs A
	void RotateMatrices(vector<SPOcd>& A, complex<double> c,
	                    complex<double> s, int i, int j);

	// Rotate Transformation matrix
	void GivensTrafoRotation(SPOcd & trafo, complex<double> c,
	                         complex<double> s, int i, int j);

	// Calculate Jacobi-Angles c, s for given elemtents i, j
	void CalculateAngles(complex<double>& c, complex<double>& s,
	                     int i, int j, const vector<SPOcd>& A);

	// Build the Matrix G, that is required to Calculate Jacobi-angles
	SPOcd BuildGMatrix(int i, int j, const vector<SPOcd>& A);

	// Weight the matrices in A with the weight matrix W
	void WeightMatrices(vector<SPOcd>& A, const SPOcd& W);

	// Build Hessian for RFO
	Matrixd RFO_BuildHessian(const Matrixd & preHessian,
	                         const Vectord & grad, double a);

	// Interpret a Vectord(2) as a complex<double>
	complex<double> InterpretComplex(const Vectord & vec);

	// Change of the Diagonals of a matrix under Givens-rotation
	Vectord RotatedDiagonals(const SPOcd & A, int p, int q,
	                         complex<double> c, complex<double> s);

	Matrixcd Rotate(const SPOcd& A,
			int p, int q, complex<double> c, complex<double> s);
};

