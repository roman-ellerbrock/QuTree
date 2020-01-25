#pragma once
#include "Tensor.h"
#include "PrimitiveBasis.h"

class SingleParticleOperator
{
public:
	SingleParticleOperator() = default;
	~SingleParticleOperator() = default;

	virtual void Apply(const PrimitiveBasis& grid, Tensorcd& hAcoeff,
		const Tensorcd& Acoeff)const = 0;
};

typedef SingleParticleOperator SPO;
