#pragma once
#include "SingleParticleOperator.h"

// This is the abstract form of a single particle operator
// implemented in functional paradigma
typedef function<void(const PrimitiveBasis&, Tensorcd&, const Tensorcd&)> SPOF;

class SingleParticleOperatorFunction
	: public SPO {
public:
	SingleParticleOperatorFunction() = default;

	SingleParticleOperatorFunction(const SPOF h)
		: h_(move(h)) {}

	~SingleParticleOperatorFunction() = default;

	void Apply(const PrimitiveBasis& grid, Tensorcd& hAcoeff,
		const Tensorcd& Acoeff) const override;

protected:
	SPOF h_;
};

typedef SingleParticleOperatorFunction SPOf;

