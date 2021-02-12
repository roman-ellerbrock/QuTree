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

	void initialize(double omega, double r0, double wfr0, double wfomega) override;

	void applyX(Tensorcd& uA, const Tensorcd & A)const override;
	void applyX2(Tensorcd& uA, const Tensorcd & A)const override;
	void applyP(Tensorcd& uA, const Tensorcd & A)const override;
	void applyKin(Tensorcd& uA, const Tensorcd & A)const override;
	void initSPF(Tensorcd & A)const override;

	const Vectord& getX()const override { return x_; }
	Vectord& getX() override { return x_; }

	void toGrid(Tensorcd& uA, const Tensorcd& A)const override;
	void fromGrid(Tensorcd& uA, const Tensorcd& A)const override;
	int oSQR()const override {return -1;}//added by TW 30.07.17
	bool hasDVR()const override {return true;}

protected:
	
	Vectord x_;
	Vectord p_;
	Matrixcd trafo_;

	double x0_, x1_, wfr0_, wfomega_;

	int dim_;
};

