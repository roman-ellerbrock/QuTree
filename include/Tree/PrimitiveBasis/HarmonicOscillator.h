#pragma once
#include "DVR.h"

class HarmonicOscillator
	: public DVR {
public:
	HarmonicOscillator() = default;

protected:
	[[nodiscard]] Tensorcd buildX(size_t dim)const override;
	[[nodiscard]] Tensorcd buildP(size_t dim)const override;
	[[nodiscard]] Tensorcd buildKin(size_t dim) const override;
};
