#pragma once
#include"Tensor.h"

class PrimitiveBasis
{
public:
	PrimitiveBasis() = default;
	virtual ~PrimitiveBasis() = default;

	virtual void Initialize(double par0, double par1, double par2, double par3) = 0;
	virtual void InitSPF(Tensorcd& Acoeffs)const = 0;
	virtual Tensorcd applyX(const Tensorcd& Acoeffs)const = 0;
	virtual Tensorcd ApplyX2(const Tensorcd& Acoeffs)const = 0;
	virtual Tensorcd ApplyP(const Tensorcd& Acoeffs)const = 0;
	virtual Tensorcd ApplyKin(const Tensorcd& Acoeffs)const = 0;
	virtual const Vectord& GetX()const = 0;
	virtual Vectord& GetX() = 0;
	virtual Tensorcd ToGrid(const Tensorcd& Acoeffs)const = 0;
	virtual Tensorcd FromGrid(const Tensorcd& Acoeffs)const = 0;
	virtual int oSQR()const = 0;
	virtual bool HasDVR()const = 0; // Tells wether a primitive basis does have a grid representation
};

