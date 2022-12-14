#pragma once
#include "TreeOperators/FortranSOP.h"

class CH3_quasiexact :
	public FortranSOP
{
public:
	CH3_quasiexact(const Tree& tree);
	~CH3_quasiexact() = default;

private:
	void callHinit(Vectorcd& coeffs, Matrix<int>& diag) override;

	void InitOperator() override;
};

