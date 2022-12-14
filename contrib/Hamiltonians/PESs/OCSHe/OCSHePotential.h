#pragma once
#include "Potential.h"

class OCSHePotential
	:public Potential
{
public:
  OCSHePotential();
	~OCSHePotential();

	void Initialize();

	double Evaluate(const Vectord& Xv, size_t part);

private:
	double cm = 219474.63068;
};

