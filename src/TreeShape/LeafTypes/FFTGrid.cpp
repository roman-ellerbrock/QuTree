#include "TreeShape/LeafTypes/FFTGrid.h"
#include "Util/FFT.h"
#include "Util/QMConstants.h"

FFTGrid::FFTGrid(int dim)
	: dim_(dim), x_(dim), p_(dim), trafo_(dim, dim) {}

void FFTGrid::initialize(double x0, double x1, double wfr0, double wfomega) {
	x0_ = x0;
	x1_ = x1;
	wfr0_ = wfr0;
	wfomega_ = wfomega;

	assert(dim_ > 0);
	// Set x_ values
	/*	double dx = (x1_ - x0_) / (dim_ - 1);
		cout << "dx= " << dx << endl;
		cout << "dim_ = " << dim_ << endl;
		for (int i = 0; i < dim_; i++) {
			x_(i) = x0_ + i * dx;
		}*/
	double dx = (x1_ - x0_) / dim_;
	for (int i = 0; i < dim_; i++) {
		x_(i) = x0 + (i + 0.5) * dx;
	}

	// Set p_ values
	/*	double pi = 3.1415926538375;
		double dp = 2 * pi / (dx * dim_);
		double prange = (dim_ - 1) * dp;
		double p0 = -prange / 2.;
		for (int i = 0; i < dim_; i++) {
			p_(i) = p0 + i * dp;
		}*/
	double dp = QM::two_pi / (dx * dim_);
	double prange = dim_ * dp;
	double p0 = -prange / 2.;
	for (int i = 0; i < dim_; i++) {
		p_(i) = p0 + (i + 0.5) * dp;
	}

	complex<double> imag(0., 1.);
	/*	for (int i = 0; i < dim_; i++)
			for (int j = 0; j < dim_; j++)
				trafo_(j, i) = exp(-imag * (x_(i) - x0_) * p_(j)) / sqrt(1. * dim_);
	//			trafo_(j, i) = exp(-imag*(x_(i) - x0_)*(p_(j) - p0));
	 */
	for (int i = 0; i < dim_; i++)
		for (int j = 0; j < dim_; j++)
			trafo_(j, i) = exp(-imag * x_(i) * p_(j)) / sqrt((double) dim_);
}
/*
void FFTGrid::initialize(double x0, double x1, double wfr0, double wfomega) {
	x0_ = x0;
	x1_ = x1;
	wfr0_ = wfr0;
	wfomega_ = wfomega;

	assert(dim_ > 0);
	// Set x_ values
	double dx = (x1_ - x0_) / (dim_ - 1);
	for (int i = 0; i < dim_; i++) {
		x_(i) = x0_ + i * dx;
	}

	// Set p_ values
	double pi = 3.1415926538375;
	double dp = 2 * pi / (dx * dim_);
	double prange = (dim_ - 1) * dp;
	double p0 = -prange / 2.;
	for (int i = 0; i < dim_; i++) {
		p_(i) = p0 + i * dp;
	}

	complex<double> imag(0., 1.);
	for (int i = 0; i < dim_; i++)
		for (int j = 0; j < dim_; j++)
			trafo_(j, i) = exp(-imag * (x_(i) - x0_) * p_(j)) / sqrt(1. * dim_);
//			trafo_(j, i) = exp(-imag*(x_(i) - x0_)*(p_(j) - p0));
}*/

void FFTGrid::applyX(Tensorcd& xA, const Tensorcd& Acoeffs) const {
//	#pragma omp for
	xA = Acoeffs;
	for (int n = 0; n < Acoeffs.shape().lastDimension(); n++)
		for (int i = 0; i < dim_; i++)
			xA(i, n) *= x_(i);
}

void FFTGrid::applyX2(Tensorcd& xA, const Tensorcd& Acoeffs) const {
//	#pragma omp for
	xA = Acoeffs;
	for (int n = 0; n < Acoeffs.shape().lastDimension(); n++)
		for (int i = 0; i < dim_; i++)
			xA(i, n) *= x_(i) * x_(i);
}

void FFTGrid::applyP(Tensorcd& pA, const Tensorcd& Acoeffs) const {
	fromGrid(pA, Acoeffs);
//	#pragma omp for
	for (int n = 0; n < Acoeffs.shape().lastDimension(); n++)
		for (int i = 0; i < dim_; i++)
			pA(i, n) *= -p_(i);

	Tensorcd tmp(pA);
	toGrid(pA, tmp);
}

void FFTGrid::toGrid(Tensorcd& uA, const Tensorcd& Acoeffs) const {
	uA = tMatrixTensor(trafo_, Acoeffs, 0);
}

void FFTGrid::fromGrid(Tensorcd& uA, const Tensorcd& Acoeffs) const {
	uA = matrixTensor(trafo_, Acoeffs, 0);
}

void FFTGrid::applyKin(Tensorcd& pA, const Tensorcd& Acoeffs) const {
	fromGrid(pA, Acoeffs);
	for (int n = 0; n < Acoeffs.shape().lastDimension(); n++)
		for (int i = 0; i < dim_; i++)
			pA(i, n) *= 0.5 * p_(i) * p_(i);
	Tensorcd tmp(pA);
	toGrid(pA, tmp);
}

void FFTGrid::initSPF(Tensorcd& phi) const {
	for (int i = 0; i < dim_; i++) {
		phi(i, 0) = exp(-0.5 * wfomega_ * pow(x_(i) - wfr0_, 2));
	}

	// excitations
	int nstates = phi.shape().lastDimension();
	for (int n = 1; n < nstates; n++) {
		for (int i = 0; i < dim_; i++) {
			phi(i, n) = x_(i) * phi(i, n - 1);
		}
		gramSchmidt(phi);
	}

	// orthonormalize
	gramSchmidt(phi);
}
