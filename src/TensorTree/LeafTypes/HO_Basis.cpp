#include "HO_Basis.h"

void HO_Basis::Initialize(double omega_, double r0_, double wfr0_, double wfomega_) {
	// save all parameters
	omega = omega_;
	r0 = r0_;
	wfr0 = wfr0_;
	wfomega = wfomega_;

	// Initialize location matrix
	FactorMatrixd xmat = InitXmat();

	Matrixcd x2(dim, dim);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			x2(j, i) = xmat(j, i);

	xmat.rDiag(trafo, x);
	for (int i = 0; i < dim; i++)
		x(i) += r0;

	// Init kinetic Energy in FBR and transform operator to DVR
	kin = InitKin();
	kin = SPOUnitarySimilarityTrafo(kin, trafo);

	// Init momentum
	p = InitPmat();
	FactorMatrixcd trafocd(trafo.Dim(), trafo.Mode());
	for (int i = 0; i < trafocd.Dim1(); i++)
		for (int j = 0; j < trafocd.Dim2(); j++)
			trafocd(j, i) = trafo(j, i);

	p = SPOUnitarySimilarityTrafo(p, trafocd);
}

void HO_Basis::InitSPF(Tensorcd& phi) const {
	TensorDim tdim(phi.Dim());
	int nstates = tdim.getntensor();
	// soft check for bottom layer_
	assert(tdim.F() == 1);
	assert(tdim.getdimpart() == dim);

	// set ground state wf
	for (int i = 0; i < dim; i++) {
//		phi(i, 0) = exp(-0.5*wfomega*pow(x(i) - wfr0, 2));
		double w = trafo(0, i) / exp(-0.5 * omega * pow(x(i) - r0, 2));
		phi(i, 0) = w * exp(-0.5 * wfomega * pow(x(i) - wfr0, 2));
	}

	// excitations
	for (int n = 1; n < nstates; n++) {
		for (int i = 0; i < dim; i++) {
			phi(i, n) = x(i) * phi(i, n - 1);
		}
	}
	// orthonormalize
	GramSchmidt(phi);
}

FactorMatrixd HO_Basis::InitXmat() {
	FactorMatrixd x(dim, 0);
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				x(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				x(m, n) = sqrt(n * 1.);
			}
			x(m, n) *= sqrt(1. / (2. * omega));
		}
	}
	return x;
}

FactorMatrixcd HO_Basis::InitPmat() {
	FactorMatrixcd p(dim, 0);
	complex<double> imag(0., 1.);
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				p(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				p(m, n) = -sqrt(n * 1.);
			}
			p(m, n) *= imag * sqrt(omega / 2.);
		}
	}
	return p;
}

FactorMatrixd HO_Basis::InitKin() {
	FactorMatrixd Kin(dim, 0);

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i == j + 2) {
				Kin(j, i) = sqrt((j + 1.) * (j + 2.));
			}
			if (i == j) {
				Kin(j, i) = -(2 * j + 1);
			}
			if (i == j - 2) {
				Kin(j, i) = sqrt(j * (j - 1.));
			}
			// factor from p-representation and from kinetic energy kin=0.5*p�
			Kin(j, i) *= (-omega / 2.) * (1 / 2.);
		}
	}
	return Kin;
}
