#pragma once
#include "Tree/PrimitiveBasis/PrimitiveBasis.h"

namespace JordanWigner {
	Tensorcd sigmaX();
	Tensorcd identity();
	Tensorcd sigmaY();
	Tensorcd sigmaZ();
	Tensorcd sigmaPlus();
	Tensorcd sigmaMinus();
}

class SpinGroup : public PrimitiveBasis {
public:
	~SpinGroup() = default;

	void initialize(const BasisParameters& par) override;
	void occupy(Tensorcd& A) const override;

	void applyX(Tensorcd& xA, const Tensorcd& A) const override { xA = A; }
	void applyX2(Tensorcd& x2A, const Tensorcd& A) const override { x2A = A; }
	void applyP(Tensorcd& pA, const Tensorcd& A) const override { pA = A; }
	void applyKin(Tensorcd& kinA, const Tensorcd& A) const override { kinA = A; }

	[[nodiscard]] const Tensord& getX() const override { cout << "Error in SpinGroup.h"; exit(1);  }
	Tensord& getX() override { cout << "Error in SpinGroup.h"; exit(1); }

	void toGrid(Tensorcd& UA, const Tensorcd& A) const override { UA = A; }
	void fromGrid(Tensorcd& UA, const Tensorcd& A) const override { UA = A; }

	[[nodiscard]] bool hasDVR() const override { return false; }

private:
	double alpha_;
};

