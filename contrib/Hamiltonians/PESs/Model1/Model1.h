#pragma once
#include "Potential.h"
class Model1 :
	public Potential
{
public:
	Model1();
	~Model1();

	vector<MultiParticleOperator> Build(const mctdhBasis& basis);

	double Evaluate(const Vectord& Xv, size_t part);

protected:
	bool simple;
};

