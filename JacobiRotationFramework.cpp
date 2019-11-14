#include "JacobiRotationFramework.h"


void JacobiRotationFramework::GivensRotation(SPOcd& B,
	complex<double> c, complex<double> s, int i, int j)
{

	// Perform a Givens-Rotation: R*A*R^H with R= {{c, s}, {-s*, c}}
	int dim = B.Dim();
	SPOcd Ccopy(dim, 2);

	for (int n = 0; n < dim; n++)
	{
		Ccopy(n, 0) = c*B(n, i) + s*B(n, j);
		Ccopy(n, 1) = c*B(n, j) - conj(s)*B(n, i);
	}

	for (int n = 0; n < dim; n++)
	{
		B(n, i) = Ccopy(n, 0);
		B(n, j) = Ccopy(n, 1);
	}

	for (int n = 0; n < dim; n++)
	{
		Ccopy(n, 0) = B(i, n);
		Ccopy(n, 1) = B(j, n);
	}

	for (int n = 0; n < dim; n++)
	{
		B(i, n) = c*Ccopy(n, 0) + conj(s)*Ccopy(n, 1);
		B(j, n) = c*Ccopy(n, 1) - s*Ccopy(n, 0);
	}
}

void JacobiRotationFramework::RotateMatrices(vector<SPOcd>& A,
	complex<double> c, complex<double> s, int i, int j)
{
	// Perform a Givens-Rotation: R*A*R^H with R= {{c, s}, {-s*, c}}
	int dim = A[0].Dim();
	SPOcd Ccopy(dim, 2);

	int nmat = A.size();
	for (int k = 0; k < nmat; k++)
	{
		SPOcd& B = A[k];

		GivensRotation(B, c, s, i, j);
	}
}

void JacobiRotationFramework::GivensTrafoRotation(SPOcd& trafo,
	complex<double> c, complex<double> s, int i, int j)
{
	int dim = trafo.Dim();
	SPOcd copy(dim, 2);

	for (int n = 0; n < dim; n++)
	{
		copy(n, 0) = trafo(n, i);
		copy(n, 1) = trafo(n, j);
	}

	for (int n = 0; n < dim; n++)
	{
		trafo(n, i) = c*copy(n, 0) + s*copy(n, 1);
		trafo(n, j) = c*copy(n, 1) - conj(s)*copy(n, 0);
	}

}

void JacobiRotationFramework::CalculateAngles(complex<double>& c,
	complex<double>& s,	int i, int j, const vector<SPOcd>& A)
{
	// Build the G-Matrix
	SPOcd G = BuildGMatrix(i, j, A);

	// Diagonalize it (phase convention is important here (x>0)!)
	SPOcd trafo(3, 3);
	Vectord  ev(3);
	G.cDiag(trafo, ev);

	// Calculate the angles from the
	complex<double> x = trafo(0, 2);
	complex<double> y = trafo(1, 2);
	complex<double> z = trafo(2, 2);

	complex<double> r = sqrt(x*x + y*y + z*z);
	complex<double> imag(0, 1);

	// Finally calculate the angles
	complex<double> t = sqrt(2.*r*(x + r));
	c = (x + r) / t;
	s = (y - imag*z) / t;
}

SPOcd JacobiRotationFramework::BuildGMatrix(int i, int j, const vector<SPOcd>& A)
{
	SPOcd G(3, 3);
	Vectorcd h(3);
	complex<double> imag(0, 1);

	for (int k = 0; k < A.size(); k++)
	{
		const SPOcd& B = A[k];

		// Build h-Vector for this Matrix
		h(0) = B(i, i) - B(j, j);
		h(1) = B(i, j) + B(j, i);
		h(2) = imag * (B(j, i) - B(i, j));

		// Add to G-Matrix
		for (int n = 0; n < 3; n++)
		{
			for (int m = 0; m < 3; m++)
			{
				G(m, n) += real(conj(h(m))*h(n));
			}
		}
	}

	return G;
}

void JacobiRotationFramework::WeightMatrices(vector<SPOcd>& A, const SPOcd & W)
{
	// Weight every Matrix in vector A with the weight matrix W
	// A = 0.5 * (W*X + X*W)
	for (int k = 0; k < A.size(); k++)
	{
		SPOcd& X = A[k];
		int mode = X.Mode();
		SPOcd Xw = 0.5*(X*W + W*X);
		A[k] = SPOcd(Xw, mode);
	}
}

Matrixd JacobiRotationFramework::RFO_BuildHessian(const Matrixd& preHessian,
	const Vectord& grad, double a)
{
	Matrixd Hessian(3, 3);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Hessian(j, i) = preHessian(j, i)*a*a;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		Hessian(2, i) = grad(i) * a;
		Hessian(i, 2) = grad(i) * a;
	}
	return Hessian;
}

Vectord JacobiRotationFramework::RotatedDiagonals(const SPOcd& A,
		int p, int q, complex<double> c, complex<double> s)
{
	// First rotation W*J^H
	Matrixcd copy(2, 2);
	copy(0, 0) = c*A(p, p) + s*A(p, q);
	copy(0, 1) = c*A(p, q) - conj(s)*A(p, p);
	copy(1, 0) = c*A(q, p) + s*A(q, q);
	copy(1, 1) = c*A(q, q) - conj(s)*A(q, p);
	//Second rotation J*(A*J^H)
	Vectord A_new(2);
	A_new(0) = real(c*copy(0, 0) + conj(s)*copy(1, 0));
	A_new(1) = real(c*copy(1, 1) - s*copy(0, 1));

	return A_new;
}

Matrixcd JacobiRotationFramework::Rotate(const SPOcd& A,
		int p, int q, complex<double> c, complex<double> s)
{
	// First rotation W*J^H
	Matrixcd copy(2, 2);
	copy(0, 0) = c*A(p, p) + s*A(p, q);
	copy(0, 1) = c*A(p, q) - conj(s)*A(p, p);
	copy(1, 0) = c*A(q, p) + s*A(q, q);
	copy(1, 1) = c*A(q, q) - conj(s)*A(q, p);
	//Second rotation J*(A*J^H)
	Matrixcd A_new(2, 2);
	A_new(0, 0) = real(c*copy(0, 0) + conj(s)*copy(1, 0));
	A_new(0, 1) = c*copy(0, 1) + conj(s)*copy(1, 1);
	A_new(1, 0) = c*copy(1, 0) - s*copy(0, 0);
	A_new(1, 1) = real(c*copy(1, 1) - s*copy(0, 1));

	return A_new;
}

complex<double> JacobiRotationFramework::InterpretComplex(const Vectord& vec)
{
	// Cast a Vectord(2) to complex<double>
	assert(vec.Dim() == 2);
	complex<double> a(vec(0), vec(1));
	return a;
}
