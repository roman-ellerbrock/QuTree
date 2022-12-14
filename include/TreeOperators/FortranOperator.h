#pragma once
#include "TreeOperators/LeafOperator.h"

typedef function<void (int*, int*, double*, double*, int*, double*, double*, double*)> FortranSystemH;

class FortranOperator :
	public LeafOperatorcd
{
	// Wrapper for Fortran-operator. Is used and managed by FortranHamiltonian. Users do not have to change this operator-class.
	// This class can be leaved unchanged, when switching systems.
	// See FortranHamiltonian for more information.
public:
	FortranOperator() = default;
	~FortranOperator() = default;

	void Initialize(int part_, int mode_, int dim_, FortranSystemH SystemH_);

	void apply(const LeafInterface& grid, Tensorcd & hAcoeff, const Tensorcd & Acoeff) const override;

protected:
	int part, mode, dim;
	FortranSystemH SystemH;
};
