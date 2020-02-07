#include "HO_Basis.h"

void HO_Basis::Initialize(double omega, double r0, double wfr0, double wfomega) {
	// save all parameters
	omega_ = omega;
	r0_ = r0;
	wfr0_ = wfr0;
	wfomega_ = wfomega;

	// Initialize location matrix
	FactorMatrixd xmat = InitXmat();

	Matrixcd x2(dim_, dim_);
	for (int i = 0; i < dim_; i++)
		for (int j = 0; j < dim_; j++)
			x2(j, i) = xmat(j, i);

	xmat.rDiag(trafo_, x_);
	for (int i = 0; i < dim_; i++)
        x_(i) += r0;

	// Init kinetic Energy in FBR and transform operator to DVR
	kin_ = InitKin();
	kin_ = SPOUnitarySimilarityTrafo(kin_, trafo_);

	// Init momentum
	p_ = InitPmat();
	FactorMatrixcd trafocd(trafo_.Dim(), trafo_.Mode());
	for (int i = 0; i < trafocd.Dim1(); i++)
		for (int j = 0; j < trafocd.Dim2(); j++)
			trafocd(j, i) = trafo_(j, i);

	p_ = SPOUnitarySimilarityTrafo(p_, trafocd);
}

void HO_Basis::InitSPF(Tensorcd& phi) const {
	TensorDim tdim(phi.Dim());
	int nstates = tdim.GetNumTensor();
	// soft check for bottom layer_
	assert(tdim.GetOrder() == 1);
	assert(tdim.GetDimPart() == dim_);

	// set ground state_ wf
	for (int i = 0; i < dim_; i++) {
//		phi(i, 0) = exp(-0.5*wfomega_*pow(x_(i) - wfr0_, 2));
		double w = trafo_(0, i) / exp(-0.5 * omega_ * pow(x_(i) - r0_, 2));
		phi(i, 0) = w * exp(-0.5 * wfomega_ * pow(x_(i) - wfr0_, 2));
	}

	// excitations
	for (int n = 1; n < nstates; n++) {
		for (int i = 0; i < dim_; i++) {
			phi(i, n) = x_(i) * phi(i, n - 1);
		}
	}
	// orthonormalize
	GramSchmidt(phi);
}

FactorMatrixd HO_Basis::InitXmat() {
	FactorMatrixd x(dim_, 0);
	for (int n = 0; n < dim_; n++) {
		for (int m = 0; m < dim_; m++) {
			if (m == n + 1) {
				x(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				x(m, n) = sqrt(n * 1.);
			}
			x(m, n) *= sqrt(1. / (2. * omega_));
		}
	}
	return x;
}

FactorMatrixcd HO_Basis::InitPmat() {
	FactorMatrixcd p(dim_, 0);
	complex<double> imag(0., 1.);
	for (int n = 0; n < dim_; n++) {
		for (int m = 0; m < dim_; m++) {
			if (m == n + 1) {
				p(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				p(m, n) = -sqrt(n * 1.);
			}
			p(m, n) *= imag * sqrt(omega_ / 2.);
		}
	}
	return p;
}

FactorMatrixd HO_Basis::InitKin() {
	FactorMatrixd Kin(dim_, 0);

	for (int i = 0; i < dim_; i++) {
		for (int j = 0; j < dim_; j++) {
			if (i == j + 2) {
				Kin(j, i) = sqrt((j + 1.) * (j + 2.));
			}
			if (i == j) {
				Kin(j, i) = -(2 * j + 1);
			}
			if (i == j - 2) {
				Kin(j, i) = sqrt(j * (j - 1.));
			}
			// factor from p_-representation and from kinetic energy kin_=0.5*p_�
			Kin(j, i) *= (-omega_ / 2.) * (1 / 2.);
		}
	}
	return Kin;
}
