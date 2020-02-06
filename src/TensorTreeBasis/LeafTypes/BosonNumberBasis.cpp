#include "TensorTreeBasis/LeafTypes/BosonNumberBasis.h"

BosonNumberBasis::BosonNumberBasis(int dim)
	:NumberBasis(dim,false)
{
}

BosonNumberBasis::~BosonNumberBasis()
{
}

// Apply primitive x_ for several single particle functions
Tensorcd BosonNumberBasis::ApplyX2(const Tensorcd& phi)const
{
	return phi;
}
