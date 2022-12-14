#pragma once
#include "Potential.h"
class CDVRModelV2 :
	public Potential
{
public:
	CDVRModelV2();
	~CDVRModelV2();

	double Evaluate(const Vectord& Xv, size_t part);

protected:
};

