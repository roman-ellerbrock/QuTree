#pragma once
#include "Tensor.h"
#include "PrimitiveBasis.h"

template <typename T>
class SingleParticleOperator
{
public:
	SingleParticleOperator() = default;
	~SingleParticleOperator() = default;

	virtual void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff)const = 0;
};

template <typename T>
using SPO = SingleParticleOperator<T>;

typedef SingleParticleOperator<complex<double>> SPOcd;
