//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "Util/RandomMatrices.h"

namespace RandomMatrices {
	Matrixcd GUE(size_t dim, mt19937& gen) {
		Matrixcd r(dim, dim);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				r(j, i) = complex<double>(dist(gen), dist(gen));
			}
		}
		return 0.5 * (r + r.Adjoint());
	}

	SpectralDecompositioncd GUE_diag(size_t dim, mt19937& gen) {
		auto A = GUE(dim, gen);
		return Diagonalize(A);
	}

	Matrixd GOE(size_t dim, mt19937& gen) {
		Matrixd r(dim, dim);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				r(j, i) = dist(gen);
			}
		}
		return 0.5 * (r + r.Adjoint());
	}

	SpectralDecompositiond GOE_diag(size_t dim, mt19937& gen) {
		auto A = GOE(dim, gen);
		return Diagonalize(A);
	}
}

