#include "Tree/PrimitiveBasis/HarmonicOscillator.h"

Tensorcd HarmonicOscillator::buildX(size_t dim) const {
	Tensorcd x({dim, dim});
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				x(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				x(m, n) = sqrt(n * 1.);
			}
			x(m, n) *= sqrt(1. / (2. * par_.omega()));
		}
	}
	return x;
}

Tensorcd HarmonicOscillator::buildP(size_t dim) const {
	Tensorcd p({dim, dim});
	complex<double> imag(0., 1.);
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				p(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				p(m, n) = -sqrt(n * 1.);
			}
			p(m, n) *= imag * sqrt(par_.omega() / 2.);
		}
	}
	return p;
}

Tensorcd HarmonicOscillator::buildKin(size_t dim) const {
	Tensorcd Kin({dim, dim});

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
			// factor from p_-representation and from kinetic energy kin_=0.5*p_�
			Kin(j, i) *= (-par_.omega() / 2.) * (1 / 2.);
		}
	}
	return Kin;
}

Tensorcd HarmonicOscillator::buildW(size_t dim) const {
	Tensorcd w({dim});
	for (int i = 0; i < dim; i++) {
		w(i) = trafo_(0, i) / exp(-0.5 * par_.omega() * pow(x_(i) - par_.r0(), 2));
	}
	return w;
}