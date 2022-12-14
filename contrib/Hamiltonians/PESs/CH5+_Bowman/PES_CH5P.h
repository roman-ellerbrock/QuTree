#pragma once
#include "Potential.h"
class PES_CH5P :
	public Potential
{
public:
	PES_CH5P();
	~PES_CH5P();

	// Evaluate the PES
	double Evaluate(const Vectord & Xv, size_t part);
};

