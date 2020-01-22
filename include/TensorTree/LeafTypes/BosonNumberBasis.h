#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "NumberBasis.h"

class BosonNumberBasis :
public NumberBasis
{
 public:
	BosonNumberBasis(int dim);
	~BosonNumberBasis();

	Tensorcd ApplyX2(const Tensorcd & phi)const;
};

