#pragma once
#include "Core/Tensor.h"
#include "TreeOperators/Hamiltonian.h"

class qbit_x
	:public LeafOperatorcd
{
	// This SPO represents the bit-operator that equals 1 if a bit
	// is at order nu and equals 0 otherwise.
	public:
		qbit_x(size_t nu_) : nu(nu_) {}
		~qbit_x() = default;

		void apply(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi)const override;

	private:
		size_t nu;
};

namespace PauliMatrices {
	void Identity(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);

	void sigma_x(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);
	void sigma_y(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);
	void sigma_z(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);

	void bit_x(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);
	void bit_x_shift(const LeafInterface& grid, Tensorcd& psi, const Tensorcd& phi);

	size_t idx(size_t i, size_t j, size_t f);
}

