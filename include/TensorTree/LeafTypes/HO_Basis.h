#pragma once
#include "SimultaneousDiagonalization.h"
#include "DVRBasis.h"

class HO_Basis
: public DVRBasis
{
public:

	HO_Basis(int dim): DVRBasis(dim){};

	void Initialize(double omega, double r0, double wfr0, double wfomega);

	void InitSPF(Tensorcd & phi)const;


protected:
	FactorMatrixd InitXmat();
	FactorMatrixcd InitPmat();
	FactorMatrixd InitKin();

	};
