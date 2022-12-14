#pragma once
#include "Potential.h"
class PES_CH4 :
	public Potential
{
public:
	PES_CH4();
	~PES_CH4();

	void Initialize();

	double Evaluate(const Vectord & Xv, size_t part);
};

