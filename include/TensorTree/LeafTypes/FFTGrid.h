#pragma once
#include "PrimitiveBasis.h"
#include "Matrix.h"
#include "Vector.h"
#include "FactorMatrix.h"
#include "Tensor.h"

class FFTGrid :
	public PrimitiveBasis
{
public:
	FFTGrid(int dim_);
	~FFTGrid() = default;

	void Initialize(double omega, double r0, double wfr0, double wfomega);

	Tensorcd applyX(const Tensorcd & Acoeffs)const;
	Tensorcd ApplyX2(const Tensorcd & Acoeffs)const;
	Tensorcd ApplyP(const Tensorcd & Acoeffs)const;
	Tensorcd ApplyKin(const Tensorcd & Acoeffs)const;
	void InitSPF(Tensorcd & Acoeffs)const;

	const Vectord& GetX()const { return x; }
	Vectord& GetX() { return x; }

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const;
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const;
	int oSQR()const{return -1;}//added by TW 30.07.17
	bool HasDVR()const {return true;}

protected:
	
	Vectord x;
	Vectord p;
	FactorMatrixcd trafo;

	double x0, x1, wfr0, wfomega;

	int dim;
};

