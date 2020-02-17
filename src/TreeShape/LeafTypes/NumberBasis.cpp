#include "NumberBasis.h"


NumberBasis::NumberBasis(int dim, bool fermion)
	: LeafInterface() {
	dim_ = dim;
	fermion_ = fermion;

	//check basis size for fermions
	if (fermion_) assert(dim_ <= 2);
}

void NumberBasis::Initialize(double occ,
	double oSQR,
	double minimalOccupation,
	double dummy3) {
	// save all parameters
	startOcc_ = (int) occ;
	osqrtree_ = (int) oSQR;
	minOcc_ = (int) minimalOccupation;

	if (minOcc_ > startOcc_) {
		cout << "Error: Minimal occupation must not be bigger than\n"
			 << "the initial occupation.\n";
		assert(0);
	}
}

void NumberBasis::InitSPF(Tensorcd& phi) const {
	TensorDim tdim(phi.Dim());
	int nstates = tdim.LastActive();

	if (fermion_ && nstates > 2) {
		cout << "Error: there can be only two fermionic states\n"
			 << "per mode.\n";
		assert(0);
	}
	if (fermion_ && minOcc_ > 2) {
		cout << "Error: For fermions the minimal Occupation must be\n"
			 << "Zero.\n";
		assert(0);
	}

	// soft check for bottom layer_
	assert(tdim.GetOrder() == 2);
	assert(tdim.LastBefore() == dim_);
	assert(startOcc_ - minOcc_ < dim_);

	// set ground state_ wf
	for (int i = 0; i < dim_; i++) {
		phi(i, 0) = 1.e-7;
	}
	phi(startOcc_ - minOcc_, 0) = 1.0;

	// excitations
	int count = 0;
	for (int n = 1; n < nstates; n++) {
		for (int i = 0; i < dim_; i++) {
			phi(i, n) = 0.0;
		}
		if (count == startOcc_ - minOcc_) count++;
		phi(count, n) = 1.0;
		count++;
	}
	// orthonormalize
	GramSchmidt(phi);
}

Tensorcd NumberBasis::ToGrid(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	TensorDim tdim = phi.Dim();
	assert(tdim.GetOrder() == 2);

	return phi;
}

Tensorcd NumberBasis::FromGrid(const Tensorcd& phi) const {
	// soft check that its really a bottom-layer_ tensor
	TensorDim tdim = phi.Dim();
	assert(tdim.GetOrder() == 2);

	return phi;
}

Tensorcd NumberBasis::ApplyKin(const Tensorcd& phi) const {
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer_ tensor
	assert(tdim.GetOrder() == 2);

	Tensorcd psi(phi.Dim());

	int nstates = tdim.LastActive();
	int active = tdim.LastBefore();

	assert(active == dim_);

	for (int n = 0; n < nstates; n++) {
		for (int i = 0; i < active; i++) {
			psi(i, n) = (1.0 * (i + minOcc_)) * phi(i, n);
		}
	}
	return psi;
}

Tensorcd NumberBasis::ApplyP(const Tensorcd& phi) const {
	const TensorDim& tdim = phi.Dim();
	Tensorcd psi(tdim, false);

	size_t prim = tdim.LastBefore();
	size_t states = tdim.LastActive();

	for (size_t n = 0; n < states; n++) {
		psi[n * prim] = 0.0;
		for (size_t i = 1; i < prim; i++) {
			psi[n * prim + i] = sqrt(1. * (i + minOcc_)) * phi[n * prim + i - 1];
		}
	}
	return psi;
}

// Apply primitive x_ for several single particle functions
Tensorcd NumberBasis::applyX(const Tensorcd& phi) const {
	const TensorDim& tdim = phi.Dim();
	Tensorcd psi(tdim, false);

	size_t prim = tdim.LastBefore();
	size_t states = tdim.LastActive();

	for (int n = 0; n < states; n++) {
		psi[n * prim + prim - 1] = 0.0;
		for (int i = 1; i < prim; i++) {
			psi[n * prim + i - 1] = sqrt(1. * (i + minOcc_)) * phi[n * prim + i];
		}
	}
	return psi;
}
