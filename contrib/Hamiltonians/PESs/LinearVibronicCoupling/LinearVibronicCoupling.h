#pragma once
#include "Potential.h"
class LinearVibronicCoupling :
	public Potential
{
public:
	LinearVibronicCoupling();
	~LinearVibronicCoupling();

	vector<MultiParticleOperator> Build(const mctdhBasis& basis);

	double Evaluate(const Vectord& Xv, size_t part);

protected:
	bool simple;
};

