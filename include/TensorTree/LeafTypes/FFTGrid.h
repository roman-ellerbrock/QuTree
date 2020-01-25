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
	explicit FFTGrid(int dim_);
	~FFTGrid() override = default;

	void Initialize(double omega, double r0, double wfr0, double wfomega) override;

	Tensorcd applyX(const Tensorcd & A)const override;
	Tensorcd ApplyX2(const Tensorcd & A)const override;
	Tensorcd ApplyP(const Tensorcd & A)const override;
	Tensorcd ApplyKin(const Tensorcd & A)const override;
	void InitSPF(Tensorcd & A)const override;

	const Vectord& GetX()const override { return x; }
	Vectord& GetX() override { return x; }

	Tensorcd ToGrid(const Tensorcd& A)const override;
	Tensorcd FromGrid(const Tensorcd& A)const override;
	int oSQR()const override {return -1;}//added by TW 30.07.17
	bool HasDVR()const override {return true;}

protected:
	
	Vectord x;
	Vectord p;
	FactorMatrixcd trafo;

	double x0, x1, wfr0, wfomega;

	int dim;
};

