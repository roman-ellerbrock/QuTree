#include "LegendrePolynomials.h"

void LegendrePolynomials::Initialize(double omega_, double r0_, double wfr0_, double wfomega_)
{
	// save all parameters
	omega = omega_;
	r0 = r0_;
	wfr0 = wfr0_;
	wfomega = wfomega_;

	//
	FactorMatrixd xmat = InitXmat();
	// @TODO: temporary test until proven!!!
	for (size_t i = 0; i < p.Dim(); i++)
		for (size_t j = 0; j < p.Dim(); j++)
			p(j, i) = xmat(j, i);
	xmat.rDiag(trafo, x);

	for (int i = 0; i < dim; i++)
	{
		x(i) = acos(x(i));
		x(i) += r0;
		x(i) = cos(x(i));
/*
		x(i) = acos(x(i));
		cout << x(i) << endl;
		*/
//		x(i) += r0;
//		x(i) = cos(x(i));
	}

	//
	kin = InitKin();
	kin = SPOUnitarySimilarityTrafo(kin, trafo);

}


void LegendrePolynomials::InitSPF(Tensorcd & phi)const
{
	// set ground state wf
	for (int i = 0; i < dim; i++)
	{
		double w = trafo(0, i);
		double xnow = acos(x(i));
//		double xnow = x(i);
		phi(i, 0) = w*exp(-0.5*wfomega*pow(xnow - wfr0, 2));
	}

	// excited state wavefunction
	for (int n=1; n< phi.Dim().GetNumTensor(); n++)
		for (int i = 0; i < dim; i++)
		{
			phi(i, n) = phi(i, n - 1)*x(i);
		}

	// orthonormalize
	GramSchmidt(phi);
}

FactorMatrixd LegendrePolynomials::InitXmat()
{
	FactorMatrixd xmat(dim, 0);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			//			if (i == j + 1) { xmat(i, j) = j*j / (4. * j*j- 1); }
			//			if (i == j - 1) { xmat(i, j) = j*j / (2. * j + 1); }
			if (i == j + 1) { xmat(i, j) = (1.*i) / sqrt(4*i*i - 1); }
			if (i == j - 1) { xmat(i, j) = (1.*j) / sqrt((2. * j + 1)*(2.*j - 1)); }
		}
	}
	return xmat;
}

FactorMatrixd LegendrePolynomials::InitKin()
{
	FactorMatrixd kin(dim, 0);
	for (int j = 0; j < dim; j++)
	{
		kin(j, j) = 0.5*j*(j + 1);
	}
	return kin;
}
