#include "Tree/PrimitiveBasis/LegendrePolynomials.h"

Tensorcd LegendrePolynomials::buildX(size_t dim) const {
	Tensorcd xmat({dim, dim});
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i == j + 1) { xmat(i, j) = (1. * i) / sqrt(4 * i * i - 1); }
			if (i == j - 1) { xmat(i, j) = (1. * j) / sqrt((2. * j + 1) * (2. * j - 1)); }
		}
	}
	return xmat;
}

Tensorcd LegendrePolynomials::buildKin(size_t dim) const {
	Tensorcd kin({dim, dim});
	for (int j = 0; j < dim; j++) {
		kin(j, j) = 0.5 * j * (j + 1);
	}
	return kin;
}

Tensorcd LegendrePolynomials::buildW(size_t dim) const {
	Tensorcd w({dim});
	for (int i = 0; i < dim; i++) {
		double x = transformX(x_(i), true);
		x -= par_.r0();
		x = transformX(x, false);
		w(i) = trafo_(0, i) / exp(-0.5 * par_.omega() * pow(x, 2));
	}
	return w;
}
