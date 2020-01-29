#pragma once
#include "PrimitiveBasis.h"

class LogicalBasis
: public PrimitiveBasis
{
public:
	LogicalBasis() = default;
	~LogicalBasis() = default;

	void Initialize(double, double, double, double);
	void InitSPF(Tensorcd& SPF)const;

	Tensorcd applyX(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyX2(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyP(const Tensorcd& Acoeffs)const;
	Tensorcd ApplyKin(const Tensorcd& Acoeffs)const;

	const Vectord& GetX()const;
	Vectord& GetX();
	Tensorcd ToGrid(const Tensorcd& Acoeffs)const;
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const;
	int oSQR()const { return 0; }
	bool HasDVR()const { return false; } 
};

