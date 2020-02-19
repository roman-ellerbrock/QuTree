#include "FermionNumberBasis.h"

FermionNumberBasis::FermionNumberBasis(int dim_)
	: NumberBasis(dim_, true) {
	assert(dim_ <= 2);
}

Tensorcd FermionNumberBasis::ApplyX2(const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	Tensorcd psi(phi.shape());

	int active = tdim.lastBefore();
	assert(active == dim_);

	for (int n = 0; n < tdim.lastDimension(); n++) {
		psi(0, n) = 1.0;
		psi(1, n) = -1.0;
	}
	return psi;
}

Tensorcd FermionNumberBasis::pauliX(const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	Tensorcd psi(phi.shape());

	int active = tdim.lastBefore();
	assert(active == dim_);

	for (int n = 0; n < tdim.lastDimension(); n++) {
		psi(1, n) = phi(0, n);
		psi(0, n) = phi(1, n);
	}
	return psi;
}

Tensorcd FermionNumberBasis::pauliY(const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	Tensorcd psi(phi.shape());

	int active = tdim.lastBefore();
	assert(active == dim_);

	complex<double> im(0.0, 1.0);
	for (int n = 0; n < tdim.lastDimension(); n++) {
		psi(1, n) = im * phi(0, n);
		psi(0, n) = -im * phi(1, n);
	}
	return psi;
}

Tensorcd FermionNumberBasis::pauliZ(const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape();
	// check that its really a bottom-layer_ tensor
	assert(tdim.order() == 2);

	Tensorcd psi(phi.shape());

	int active = tdim.lastBefore();
	assert(active == dim_);

	for (int n = 0; n < tdim.lastDimension(); n++) {
		psi(1, n) = phi(1, n);
		psi(0, n) = -phi(0, n);
	}
	return psi;
}
