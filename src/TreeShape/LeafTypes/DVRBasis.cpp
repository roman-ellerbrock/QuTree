#include "TreeShape/LeafTypes/DVRBasis.h"

DVRBasis::DVRBasis(int dim)
	: dim_(dim), trafo_(dim, dim), x_(dim), kin_(dim, dim), p_(dim, dim),
	  omega_(0), r0_(0), wfr0_(0), wfomega_(0) {
}

void DVRBasis::applyX(Tensorcd& xPhi, const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	int active = tdim.lastBefore();
	assert(active == dim_);

	//	psi = kin_*phi;
	// @TODO: rewrite this code as a matrix*Tensor routine
//  	#pragma omp for
	for (int n = 0; n < tdim.lastDimension(); n++) {
		for (int i = 0; i < active; i++) {
			xPhi(i, n) = x_(i) * phi(i, n);
		}
	}
}

void DVRBasis::applyX2(Tensorcd& psi, const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	int active = tdim.lastBefore();
	assert(active == dim_);

	//	psi = kin_*phi;
	// @TODO: rewrite this code as a matrix*Tensor routine
//    #pragma omp for
	for (int n = 0; n < tdim.lastDimension(); n++) {
		for (int i = 0; i < active; i++) {
			psi(i, n) = x_(i) * x_(i) * phi(i, n);
		}
	}
}

void DVRBasis::applyP(Tensorcd& uPhi, const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorShape& tdim = phi.shape();
	assert(tdim.order() == 2);

	uPhi = matrixTensor(p_, phi, 0);
}

void DVRBasis::applyKin(Tensorcd& uPhi, const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorShape& tdim = phi.shape();
	assert(tdim.order() == 2);

	uPhi = matrixTensor(kin_, phi, 0);
}

void DVRBasis::ToGrid(Tensorcd& uPhi, const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorShape& tdim = phi.shape();
	assert(tdim.order() == 2);

	uPhi = multATB(trafo_, phi, 0);
}

void DVRBasis::FromGrid(Tensorcd& uPhi, const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	const TensorShape& tdim = phi.shape();
	assert(tdim.order() == 2);

	uPhi = matrixTensor(trafo_, phi, 0);
}
