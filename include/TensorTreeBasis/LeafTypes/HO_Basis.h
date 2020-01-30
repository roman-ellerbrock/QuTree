#pragma once
#include "Core/SimultaneousDiagonalization.h"
#include "DVRBasis.h"

class HO_Basis
	: public DVRBasis {
public:
	HO_Basis() = default;
	~HO_Basis() override = default;

	explicit HO_Basis(int dim)
		: DVRBasis(dim) {};

	void Initialize(double omega, double r0, double wfr0, double wfomega) override;

	void InitSPF(Tensorcd& phi) const override;


protected:
	FactorMatrixd InitXmat();
	FactorMatrixcd InitPmat();
	FactorMatrixd InitKin();
};