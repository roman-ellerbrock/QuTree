#pragma once
#include "Tree/PrimitiveBasis/DVR.h"

class FFTGrid : public DVR {
public:
	FFTGrid() = default;
	~FFTGrid() = default;

	void initialize(size_t dim, const BasisParameters& par) override;

	[[nodiscard]] Tensord buildXvec(size_t dim)const;
	[[nodiscard]] Tensorcd buildU(size_t dim)const;

	[[nodiscard]] Tensorcd buildP(size_t dim)const override;
	[[nodiscard]] Tensorcd buildKin(size_t dim) const override;
	[[nodiscard]] Tensorcd buildW(size_t dim) const override;

};

