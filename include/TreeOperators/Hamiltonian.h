#pragma once
#include "TreeOperators/SumOfProductsOperator.h"

class Hamiltonian : public SOPcd
{
public:
	Hamiltonian() : hasV(false) {}
	~Hamiltonian() = default;

	PotentialOperator V_;
	bool hasV;

	using SOPcd::operator=;

};

