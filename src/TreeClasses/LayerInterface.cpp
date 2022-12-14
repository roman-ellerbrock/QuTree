//
// Created by Roman Ellerbrock on 3/8/20.
//
#include "TreeClasses/LayerInterface.h"

double LayerInterface::Error(const Tensorcd& Phi, const Tensorcd& Chi) const {
	double Delta = 0;
	double norm = 1.;
	Tensorcd diff = Phi - Chi;

	if (node_->isToplayer()) {
		for (size_t i = 0; i < Phi.shape().totalDimension(); ++i) {
			Delta += real(conj(diff(i)) * diff(i));
		}

		auto S = Phi.dotProduct(Phi);
		norm = abs(S.trace());

	} else {
		const Matrixcd& rho = hRep_->rho_[*node_];
		for (size_t n = 0; n < Phi.shape().lastDimension(); ++n) {
			for (size_t m = 0; m < Phi.shape().lastDimension(); ++m) {
				for (size_t i = 0; i < Phi.shape().lastBefore(); ++i) {
					Delta += real(conj(diff(i, n)) * rho(n, m) * diff(i, m));
				}
			}
		}
		norm = abs(rho.trace());
	}
	return sqrt(abs(Delta / norm));
}

