#pragma once
#include "TreeOperators/CoordinateTransformation.h"


class TrafoCH3Quasiexact:
	public CoordinateTransformation {
public:
	TrafoCH3Quasiexact(Vectord mass_);

	Vectord transform(const Vectord& q) const override;

private:
	Vectord mass;
};


