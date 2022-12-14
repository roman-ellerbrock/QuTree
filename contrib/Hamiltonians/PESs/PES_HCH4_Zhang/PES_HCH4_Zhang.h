#pragma once
#include "Potential.h"
class PES_HCH4_Zhang :
	public Potential
{
public:
	PES_HCH4_Zhang();
	~PES_HCH4_Zhang();

	double Evaluate(const Vectord& Xv, size_t part);

protected:
	void Initialize();

};

