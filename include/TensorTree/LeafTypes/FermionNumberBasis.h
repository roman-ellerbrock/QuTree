#pragma once
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "NumberBasis.h"

class FermionNumberBasis :
public NumberBasis
{
 public:
	FermionNumberBasis(int dim);
	~FermionNumberBasis();

	Tensorcd ApplyX2(const Tensorcd & phi)const;
        
	Tensorcd pauliX(const Tensorcd & phi)const;
	Tensorcd pauliY(const Tensorcd & phi)const;
	Tensorcd pauliZ(const Tensorcd & phi)const;
};

