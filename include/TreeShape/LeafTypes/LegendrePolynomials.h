#pragma once
#include <iostream>
#include "DVRBasis.h"

class LegendrePolynomials
	: public DVRBasis {
public:
	explicit LegendrePolynomials(int dim)
		: DVRBasis(dim) {};
	void Initialize(double omega, double r0, double wfr0, double wfomega) override;
	void InitSPF(Tensorcd& Acoeffs) const override;

	Tensorcd ApplyP(const Tensorcd& Acoeffs) const override {
		cerr << "ApplyP of LegendrePolynomials has been called.\n";
		cerr << "This function is not defined for LegendrePolynomials.\n";
		exit(1);
	}

protected:
	Matrixd InitXmat();
	Matrixd InitKin();
};
