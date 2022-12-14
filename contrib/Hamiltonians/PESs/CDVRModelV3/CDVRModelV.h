#pragma once
#include "Potential.h"
class CDVRModelV3 :
	public Potential
{
public:
	CDVRModelV3();
	~CDVRModelV3();

	vector<MultiParticleOperator> Build(const mctdhBasis& basis);

	double Evaluate(const Vectord& Xv, size_t part);

protected:
	bool simple;
};

