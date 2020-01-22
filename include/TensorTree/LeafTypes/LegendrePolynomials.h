#pragma once
#include <iostream>
#include "DVRBasis.h"

class LegendrePolynomials
: public DVRBasis
{
public:
	LegendrePolynomials(int dim): DVRBasis(dim){};
	void Initialize(double omega, double r0, double wfr0, double wfomega);
	void InitSPF(Tensorcd & Acoeffs)const;

	Tensorcd ApplyP(const Tensorcd & Acoeffs)const
	{
		cerr << "ApplyP of LegendrePolynomials has been called.\n";
		cerr << "This function is not defined for LegendrePolynomials.\n";
		exit(1);
	}



protected:

	FactorMatrixd InitXmat();
	FactorMatrixd InitKin();
	//SPOd InitPmat()
	//{
	//	cerr << "InitPmat of LegendrePolynomials has been called.\n";
	//	cerr << "This function is not defined for LegendrePolynomials.\n";
	//	exit(1);
	//}

};
