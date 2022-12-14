#pragma once
#include "Potential.h"

class HeHePotential :
	public Potential
{
public:
	HeHePotential();
	~HeHePotential();

	double Evaluate(const double x);
};

