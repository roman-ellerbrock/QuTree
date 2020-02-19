#pragma once
#include "LeafInterface.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "Core/Tensor.h"

class FFTGrid :
	public LeafInterface
{
public:
	explicit FFTGrid(int dim);
	~FFTGrid() override = default;

	void Initialize(double omega, double r0, double wfr0, double wfomega) override;

	Tensorcd applyX(const Tensorcd & A)const override;
	Tensorcd ApplyX2(const Tensorcd & A)const override;
	Tensorcd ApplyP(const Tensorcd & A)const override;
	Tensorcd ApplyKin(const Tensorcd & A)const override;
	void InitSPF(Tensorcd & A)const override;

	const Vectord& GetX()const override { return x_; }
	Vectord& GetX() override { return x_; }

	Tensorcd ToGrid(const Tensorcd& A)const override;
	Tensorcd FromGrid(const Tensorcd& A)const override;
	int oSQR()const override {return -1;}//added by TW 30.07.17
	bool HasDVR()const override {return true;}

protected:
	
	Vectord x_;
	Vectord p_;
	Matrixcd trafo_;

	double x0_, x1_, wfr0_, wfomega_;

	int dim_;
};

