#include "Tree/PrimitiveBasis//FFTGrid.h"
#include "Tensor/Tensor"
#include "Util/QMConstants.h"

void FFTGrid::initialize(size_t dim, const BasisParameters& par) {
	par_ = par;

	/// The order is important.
	x_ = buildXvec(dim);
	p_ = buildP(dim);
	kin_ = buildKin(dim);
	trafo_ = buildU(dim);
	unitarySimilarityTrafo(p_, trafo_);
	unitarySimilarityTrafo(kin_, trafo_);
}

Tensorcd FFTGrid::buildU(size_t dim) const {
	Tensorcd U({dim, dim});
	double x0 = par_.par0_;
	complex<double> imag(0., 1.);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			U(j, i) = exp(-imag * (x_(i) - x0) * p_(j, j)) / sqrt(1. * dim);

	return U;
}

Tensord FFTGrid::buildXvec(size_t dim) const {
	Tensord x({dim});
	double x1 = par_.par1_;
	double x0 = par_.par0_;
	double dx = (x1 - x0) / (dim - 1);
	for (int i = 0; i < dim; i++) {
		x(i) = x0 + i * dx;
	}
	return x;
}

Tensorcd FFTGrid::buildP(size_t dim) const {
	Tensorcd p({dim, dim});
	// Set p_ values
	double x1 = par_.par1_;
	double x0 = par_.par0_;
	double dx = (x1 - x0) / (dim - 1);
	double dp = QM::two_pi / (dx * dim);
	double prange = (dim - 1) * dp;
	double p0 = -prange / 2.;
	for (int i = 0; i < dim; i++) {
		p(i, i) = p0 + i * dp;
	}
	return p;
}

Tensorcd FFTGrid::buildKin(size_t dim) const {
	auto T = buildP(dim);
	T = 0.5 * gemm(T, T);
	return T;
}

Tensorcd FFTGrid::buildW(size_t dim) const {
	Tensorcd w({dim});
	for (int i = 0; i < dim; i++) {
		w(i) = 1.;
	}
	return w;
}

