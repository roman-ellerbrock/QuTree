#pragma once
#include "Tensor.h"
#include "PrimitiveBasis.h"
#include "Matrix.h"
#include "Vector.h"
#include "FactorMatrix.h"

class DVRBasis
: public PrimitiveBasis
{
public:
   DVRBasis(int dim);
	~DVRBasis() = default;

	virtual void Initialize(double par0, double par1, double par2, double par3) = 0;

	virtual void InitSPF(Tensorcd& Acoeffs)const = 0;
	Tensorcd applyX(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyX2(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyP(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyKin(const Tensorcd& Acoeffs)const;

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const;
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const;
	int oSQR()const {return -1;}
	
	Vectord& GetX() {return x;}
	const Vectord& GetX() const {return x;}
	bool HasDVR()const {return true;}


protected:

	// location operators
	Vectord x;
	FactorMatrixd trafo;

	// derivative operators
//	Matrixcd p; // @TODO: implement p
	FactorMatrixd kin;
	FactorMatrixcd p;

	// parameters that define the basis
	double omega, r0, wfr0, wfomega;
	int dim;

};
