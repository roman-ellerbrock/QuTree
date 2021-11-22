#pragma once
#include <iostream>
#include "DVR.h"

class LegendrePolynomials : public DVR {
public:
	LegendrePolynomials() = default;
	~LegendrePolynomials() = default;

	[[nodiscard]] double transformX(double x, bool forth) const override {
		if (forth) {
			return acos(x);
		} else {
			return cos(x);
		}
	}


	[[nodiscard]] Tensorcd buildX(size_t dim)const override;
	[[nodiscard]] Tensorcd buildKin(size_t dim) const override;
	[[nodiscard]] Tensorcd buildW(size_t dim) const override;

};
