#pragma once
#include "PrimitiveBasis.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "Core/Tensor.h"
#include <random>

class SpinGroup :
	public PrimitiveBasis
{
public:
	explicit SpinGroup(size_t dim): dim_(dim), alpha_(0.), last_(false) {}
	~SpinGroup() = default;

	void Initialize(double par0, double par1, double par2, double par3) override;

	Tensorcd applyX(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyX2(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyP(const Tensorcd & Acoeffs)const override;
	Tensorcd ApplyKin(const Tensorcd & Acoeffs)const override;
	void InitSPF(Tensorcd & Acoeffs)const override;

	const Vectord& GetX()const override { return x_; } // does nothing
	Vectord& GetX() override{ return x_; } // does nothing

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const override { return Acoeffs; } // does nothing
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const override { return Acoeffs; } // does nothing

	int oSQR()const override {return -1;}
	bool HasDVR()const override {return false;}


private:
	Vectord x_;

	size_t dim_;
	bool last_;
	double alpha_;
};

