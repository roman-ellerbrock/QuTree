#pragma once
#include "Hamiltonian.h"

namespace scatteringoperators {
	Tensorcd Applyh(const PrimitiveBasis& grid, const Tensorcd& Psi);
	void ApplyFlux(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi);

	class CAP: public StandardRefSPO {
	public:
		CAP(double x0_, double length_, double strength_)
			: x0(x0_), length(length_), strength(strength_) {
		}

		~CAP() = default;

		void Apply(Tensorcd& HPsi, const PrimitiveBasis& grid, const Tensorcd& Psi)const override;

	protected:
		double x0;
		double length;
		double strength;
	};
}
class scattering:
	public SOP {
public:
	//! Constructor
	scattering(const mctdhBasis& basis) {
		Initialize(basis);
	}

	//! Destructor
	~scattering() = default;

	//! Initialize the Operator
	void SpecialInitialize(const mctdhBasis& basis);
};

