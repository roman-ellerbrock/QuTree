#pragma once
#include "Potential.h"
class PES_CH4SM :
	public Potential
{
public:
	PES_CH4SM();
	~PES_CH4SM();

	void Initialize();

	double Evaluate(const Vectord & Xv, size_t part);
};

