#include "DVRBasis.h"

DVRBasis::DVRBasis(int dim)
	: dim_(dim), trafo_(dim, dim), x_(dim), kin_(dim, dim), p_(dim, dim),
	  omega_(0), r0_(0), wfr0_(0), wfomega_(0) {
}

Tensorcd DVRBasis::applyX(const Tensorcd& phi) const {
	const TensorDim& tdim = phi.Dim();
	// check that its really a bottom-layer_ tensor
	assert(tdim.GetOrder() == 1);

	Tensorcd psi(phi.Dim());

	int active = tdim.GetDimPart();
	assert(active == dim_);

	//	psi = kin_*phi;
	// @TODO: rewrite this code as a matrix*Tensor routine
//  	#pragma omp for
	for (int n = 0; n < tdim.GetNumTensor(); n++) {
		for (int i = 0; i < active; i++) {
			psi(i, n) = x_(i) * phi(i, n);
		}
	}
	return psi;
}

Tensorcd DVRBasis::ApplyX2(const Tensorcd& phi) const {
	const TensorDim& tdim = phi.Dim();
	// check that its really a bottom-layer_ tensor
	assert(tdim.GetOrder() == 1);

	Tensorcd psi(phi.Dim());

	int active = tdim.GetDimPart();
	assert(active == dim_);

	//	psi = kin_*phi;
	// @TODO: rewrite this code as a matrix*Tensor routine
//    #pragma omp for
	for (int n = 0; n < tdim.GetNumTensor(); n++) {
		for (int i = 0; i < active; i++) {
			psi(i, n) = x_(i) * x_(i) * phi(i, n);
		}
	}
	return psi;
}

Tensorcd DVRBasis::ApplyP(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorDim& tdim = phi.Dim();
	assert(tdim.GetOrder() == 1);

	return multAB(p_, phi, 0);
}

Tensorcd DVRBasis::ApplyKin(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorDim& tdim = phi.Dim();
	assert(tdim.GetOrder() == 1);

	return multAB(kin_, phi, 0);
}

Tensorcd DVRBasis::ToGrid(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorDim& tdim = phi.Dim();
	assert(tdim.GetOrder() == 1);

	return multATB(trafo_, phi, 0);
}

Tensorcd DVRBasis::FromGrid(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorDim& tdim = phi.Dim();
	assert(tdim.GetOrder() == 1);

	return multAB(trafo_, phi, 0);
}
