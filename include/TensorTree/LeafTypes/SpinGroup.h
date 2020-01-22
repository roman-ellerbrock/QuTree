#pragma once
#include "PrimitiveBasis.h"
#include "Matrix.h"
#include "Vector.h"
#include "Tensor.h"
#include <random>

class SpinGroup :
	public PrimitiveBasis
{
public:
	explicit SpinGroup(size_t dim_):dim(dim_), alpha(0.), last(false) {}
	~SpinGroup() = default;

	void Initialize(double par0, double par1, double par2, double par3) override;

	Tensorcd applyX(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyX2(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyP(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyKin(const Tensorcd & Acoeffs)const override;
	void InitSPF(Tensorcd & Acoeffs)const override;

	const Vectord& GetX()const override { return x; } // does nothing
	Vectord& GetX() override{ return x; } // does nothing

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const override { return Acoeffs; } // does nothing
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const override { return Acoeffs; } // does nothing

	int oSQR()const override {return -1;}
	bool HasDVR()const override {return false;}


private:
	Vectord x;

	size_t dim;
	bool last;
	double alpha;
};

