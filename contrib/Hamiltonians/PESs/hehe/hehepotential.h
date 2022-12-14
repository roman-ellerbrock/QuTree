#pragma once
#include "stdafx.h"

class HeHePotential
{
public:
	HeHePotential();
	~HeHePotential();

	double Evaluate(const double x, size_t part);
};

