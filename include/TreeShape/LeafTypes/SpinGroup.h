#pragma once
#include "LeafInterface.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "Core/Tensor.h"
#include <random>

class SpinGroup:
	public LeafInterface {
public:
	explicit SpinGroup(size_t dim)
		: dim_(dim), alpha_(0.), last_(false) {}

	~SpinGroup() = default;

	void initialize(double par0, double par1, double par2, double par3) override;

	void applyX(Tensorcd& uA, const Tensorcd& A) const override { uA = A; }

	void applyX2(Tensorcd& uA, const Tensorcd& A) const override { uA = A; }

	void applyP(Tensorcd& uA, const Tensorcd& A) const override { uA = A; }

	void applyKin(Tensorcd& uA, const Tensorcd& A) const override { uA = A; }

	void initSPF(Tensorcd& Acoeffs) const override;

	const Vectord& getX() const override { return x_; } // does nothing
	Vectord& getX() override { return x_; } // does nothing

	void toGrid(Tensorcd& uA, const Tensorcd& Acoeffs) const override { uA = Acoeffs; } // does nothing
	void fromGrid(Tensorcd& uA, const Tensorcd& Acoeffs) const override { uA = Acoeffs; } // does nothing

	int oSQR() const override { return -1; }

	bool hasDVR() const override { return false; }


private:
	Vectord x_;

	size_t dim_;
	bool last_;
	double alpha_;
};

