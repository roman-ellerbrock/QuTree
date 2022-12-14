#pragma once
#include "Tensor.h"
#include "Hamiltonian.h"

class qbit_x
	:public BottomLayerSPO
{
	// This SPO represents the bit-operator that equals 1 if a bit
	// is at order nu and equals 0 otherwise.
	public:
		qbit_x(size_t nu_) : nu(nu_) {}
		~qbit_x() = default;

		void Apply(Tensorcd& xPhi, const PrimitiveBasis& grid, const Tensorcd& phi)const override;

	private:
		size_t nu;
};

namespace PauliMatrices {
	void Identity(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);

	void sigma_x(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);
	void sigma_y(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);
	void sigma_z(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);

	void bit_x(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);
	void bit_x_shift(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi);

	size_t idx(size_t i, size_t j, size_t f);
}

