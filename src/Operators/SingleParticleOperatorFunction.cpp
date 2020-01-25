#include "SingleParticleOperatorFunction.h"


void SPOf::Apply(const PrimitiveBasis& grid,
	Tensorcd& hAcoeff, const Tensorcd& Acoeff) const {
	h_(hAcoeff, grid, Acoeff);
}

