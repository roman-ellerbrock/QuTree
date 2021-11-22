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

	void initialize(size_t dim, const BasisParameters& par) override;
	void occupy(Tensorcd& A) const override;

private:
	double alpha_;
};

