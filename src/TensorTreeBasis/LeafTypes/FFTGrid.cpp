#include "FFTGrid.h"

FFTGrid::FFTGrid(int dim_)
	:dim(dim_), x(dim_), p(dim_), trafo(dim_, 0) { }

void FFTGrid::Initialize(double x0_, double x1_, double wfr0_, double wfomega_)
{
	x0 = x0_;
	x1 = x1_;
	wfr0 = wfr0_;
	wfomega = wfomega_;

	assert(dim > 0);
	// Set x values
	double dx = (x1 - x0) / (dim - 1);
	cout << "dx= " << dx << endl;
	cout << "dim = " << dim << endl;
	for (int i = 0; i < dim; i++)
	{
		x(i) = x0 + i*dx;
	}

	// Set p values
	double pi = 3.1415926538375;
	double dp = 2 * pi / (dx*dim);
	double prange = (dim - 1)*dp;
	double p0 = -prange / 2.;
	for (int i = 0; i < dim; i++)
	{
		p(i) = p0 + i*dp;
	}

	complex<double> imag(0., 1.);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			trafo(j, i) = exp(-imag*(x(i) - x0)*p(j)) / sqrt(1.*dim);
//			trafo(j, i) = exp(-imag*(x(i) - x0)*(p(j) - p0));
}

Tensorcd FFTGrid::applyX(const Tensorcd& Acoeffs)const
{
	Tensorcd xA(Acoeffs);
//	#pragma omp for
	for (int n = 0; n < Acoeffs.Dim().GetNumTensor(); n++)
		for (int i = 0; i < dim; i++)
			xA(i, n) *= x(i);
	return xA;
}

Tensorcd FFTGrid::ApplyX2(const Tensorcd& Acoeffs)const
{
	Tensorcd xA(Acoeffs);
//	#pragma omp for
	for (int n = 0; n < Acoeffs.Dim().GetNumTensor(); n++)
		for (int i = 0; i < dim; i++)
			xA(i, n) *= x(i)*x(i);
	return xA;
}

Tensorcd FFTGrid::ApplyP(const Tensorcd& Acoeffs)const
{
	Tensorcd pA = FromGrid(Acoeffs);
//	#pragma omp for
	for (int n = 0; n < Acoeffs.Dim().GetNumTensor(); n++)
		for (int i = 0; i < dim; i++)
			pA(i, n) *= -p(i);
	pA = ToGrid(pA);
	return pA;
}

#include "Core/FFT.h"
Tensorcd FFTGrid::ToGrid(const Tensorcd& Acoeffs)const {
	return multATB(trafo, Acoeffs);
}

Tensorcd FFTGrid::FromGrid(const Tensorcd& Acoeffs)const {
	return trafo*Acoeffs;
}

Tensorcd FFTGrid::ApplyKin(const Tensorcd& Acoeffs)const {
    Tensorcd pA = FromGrid(Acoeffs);
	for (int n = 0; n < Acoeffs.Dim().GetNumTensor(); n++)
		for (int i = 0; i < dim; i++)
			pA(i, n) *= 0.5*p(i)*p(i);
    return ToGrid(pA);
}

void FFTGrid::InitSPF(Tensorcd& phi)const {
	for (int i = 0; i < dim; i++) {
		phi(i, 0) = exp(-0.5*wfomega*pow(x(i) - wfr0, 2));
	}

	// excitations
	int nstates = phi.Dim().GetNumTensor();
	for (int n = 1; n < nstates; n++) {
		for (int i = 0; i < dim; i++) {
			phi(i, n) = x(i)*phi(i, n - 1);
		}
	}

	// orthonormalize
	GramSchmidt(phi);
}
