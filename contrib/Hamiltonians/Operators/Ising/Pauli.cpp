#include "Pauli.h"

void qbit_x::Apply(Tensorcd& xPhi, const PrimitiveBasis& grid, const Tensorcd& Phi) const {
	const TensorDim& tdim = Phi.Dim();
	size_t ntensor = tdim.getntensor();
	size_t dimpart = tdim.getdimpart();
	xPhi.Zero();
	for (size_t n = 0; n < ntensor; ++n) {
		xPhi(nu, n) = Phi(nu, n);
	}
}

namespace PauliMatrices {
	void Identity(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(0, n) = phi(0, n);
			psi(1, n) = phi(1, n);
		}
	}

	void sigma_x(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(1, n) = phi(0, n);
			psi(0, n) = phi(1, n);
		}
	}

	void sigma_y(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		complex<double> im(0.0, 1.0);
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(1, n) = im * phi(0, n);
			psi(0, n) = -im * phi(1, n);
		}
	}

	void sigma_z(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(1, n) = phi(1, n);
			psi(0, n) = -phi(0, n);
		}
	}

	void bit_x(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(1, n) = phi(1, n);
			psi(0, n) = 0;
		}
	}

	void bit_x_shift(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi) {
		const TensorDim& tdim = phi.Dim();
		constexpr size_t N = 5;
		for (int n = 0; n < tdim.getntensor(); n++) {
			psi(1, n) = -phi(1, n) + 1. / ((double) N);
			psi(0, n) = 1. / ((double) N);
		}
	}

	void bit_xs(Tensorcd& xPhi, const PrimitiveBasis& grid, const Tensorcd& phi) {
	}

	size_t idx(size_t i, size_t j, size_t f) {
		// map 2D indices (i, j) to superindex I (periodic)
		size_t red_i = i % f;
		size_t red_j = j % f;
		return red_j * f + red_i;
	}
}

