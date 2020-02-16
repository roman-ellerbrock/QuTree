#include "BosonNumberBasis.h"

BosonNumberBasis::BosonNumberBasis(int dim_)
	: NumberBasis(dim_, false) {
}

// Apply primitive x_ for several single particle functions
Tensorcd BosonNumberBasis::ApplyX2(const Tensorcd& phi) const {
	return phi;
}
