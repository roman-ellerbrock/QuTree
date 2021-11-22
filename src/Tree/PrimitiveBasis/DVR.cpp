#include "Tree/PrimitiveBasis/DVR.h"
#include "Tensor/Tensor"

void DVR::initialize(size_t dim, const BasisParameters& par) {
	par_ = par;

	Tensorcd xmat = buildX(dim);
	auto diag = heev(xmat);
	trafo_ = diag.U();
	x_ = diag.ev();

	shift(x_, par.r0());

	kin_ = buildKin(dim);
	kin_ = unitarySimilarityTrafo(kin_, trafo_);

	p_ = buildP(dim);
	p_ = unitarySimilarityTrafo(p_, trafo_);

	w_ = buildW(dim);
}

void DVR::shift(Tensord& x, double delta) const {
	for (size_t i = 0; i < x.shape_.totalDimension(); ++i) {
		x(i) = transformX(x(i), true);
		x(i) += delta;
		x(i) = transformX(x(i), false);
	}
}


void DVR::occupy(Tensorcd& phi) const {
	// set ground state wf
	const TensorShape& shape = phi.shape_;
	for (int i = 0; i < shape.lastBefore(); i++) {
		double x = transformX(x_(i), true) - par_.wfr0();
		phi(i, 0) = w_(i) * exp(-0.5 * par_.wfomega() * pow(x, 2));
	}

	// excitations
	for (int n = 1; n < shape.lastDimension(); n++) {
		for (int i = 0; i < shape.lastBefore(); i++) {
			phi(i, n) = x_(i) * phi(i, n - 1);
		}
	}

	// orthonormalize
	gramSchmidt(phi);
}

void DVR::applyX(Tensorcd& xPhi, const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape_;
	assert(tdim.order() == 2);
	assert(x_.shape_[0] == tdim.lastBefore());

	for (size_t n = 0; n < tdim.lastDimension(); n++) {
		for (size_t i = 0; i < tdim.lastBefore(); i++) {
			xPhi(i, n) = x_(i) * phi(i, n);
		}
	}
}

void DVR::applyX2(Tensorcd& psi, const Tensorcd& phi) const {
	const TensorShape& tdim = phi.shape_;
	assert(tdim.order() == 2);
	assert(x_.shape_.order() == 2);
	assert(x_.shape_[0] == tdim.lastBefore());

	for (size_t n = 0; n < tdim.lastDimension(); n++) {
		for (size_t i = 0; i < tdim.lastBefore(); i++) {
			psi(i, n) = x_(i) * x_(i) * phi(i, n);
		}
	}
}

void DVR::applyP(Tensorcd& uPhi, const Tensorcd& phi) const {
	uPhi = matrixTensor(p_, phi, 0);
}

void DVR::applyKin(Tensorcd& uPhi, const Tensorcd& phi) const {
	uPhi = matrixTensor(kin_, phi, 0);
}

void DVR::toGrid(Tensorcd& uPhi, const Tensorcd& phi) const {
	uPhi = matrixTensor(adjoint(trafo_), phi, 0);
}

void DVR::fromGrid(Tensorcd& uPhi, const Tensorcd& phi) const {
	uPhi = matrixTensor(trafo_, phi, 0);
}

Tensorcd DVR::buildW(size_t dim) const {
	Tensorcd w({dim});
	for (int i = 0; i < dim; i++) {
		double x = transformX(x_(i), true);
		x -= par_.r0();
		x = transformX(x, false);
		w(i) = trafo_(0, i) / exp(-0.5 * par_.omega() * pow(x, 2));
	}
	return w;
}