#include "TreeShape/LeafTypes/LegendrePolynomials.h"

void LegendrePolynomials::Initialize(double omega, double r0, double wfr0, double wfomega) {
	// save all parameters
	omega_ = omega;
	r0_ = r0;
	wfr0_ = wfr0;
	wfomega_ = wfomega;

	//
	Matrixd xmat = InitXmat();
	// @TODO: temporary test until proven!!!
	for (size_t i = 0; i < p_.Dim2(); i++)
		for (size_t j = 0; j < p_.Dim1(); j++)
			p_(j, i) = xmat(j, i);
	xmat.rDiag(trafo_, x_);

	for (int i = 0; i < dim_; i++) {
		x_(i) = acos(x_(i));
		x_(i) += r0;
		x_(i) = cos(x_(i));
/*
		x_(i) = acos(x_(i));
		cout << x_(i) << endl;
		*/
//		x_(i) += r0_;
//		x_(i) = cos(x_(i));
	}

	//
	kin_ = InitKin();
	kin_ = UnitarySimilarityTrafo(kin_, trafo_);
}

void LegendrePolynomials::InitSPF(Tensorcd& phi) const {
	// set ground state_ wf
	for (int i = 0; i < dim_; i++) {
		double w = trafo_(0, i);
		double xnow = acos(x_(i));
//		double xnow = x_(i);
		phi(i, 0) = w * exp(-0.5 * wfomega_ * pow(xnow - wfr0_, 2));
	}

	// excited state_ wavefunction
	for (int n = 1; n < phi.shape().lastDimension(); n++)
		for (int i = 0; i < dim_; i++) {
			phi(i, n) = phi(i, n - 1) * x_(i);
		}

	// orthonormalize
	GramSchmidt(phi);
}

Matrixd LegendrePolynomials::InitXmat() {
	Matrixd xmat(dim_, dim_);
	for (int i = 0; i < dim_; i++) {
		for (int j = 0; j < dim_; j++) {
			//			if (i == j + 1) { xmat(i, j) = j*j / (4. * j*j- 1); }
			//			if (i == j - 1) { xmat(i, j) = j*j / (2. * j + 1); }
			if (i == j + 1) { xmat(i, j) = (1. * i) / sqrt(4 * i * i - 1); }
			if (i == j - 1) { xmat(i, j) = (1. * j) / sqrt((2. * j + 1) * (2. * j - 1)); }
		}
	}
	return xmat;
}

Matrixd LegendrePolynomials::InitKin() {
	Matrixd kin(dim_, dim_);
	for (int j = 0; j < dim_; j++) {
		kin(j, j) = 0.5 * j * (j + 1);
	}
	return kin;
}
