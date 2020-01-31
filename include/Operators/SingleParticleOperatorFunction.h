#pragma once
#include "SingleParticleOperator.h"

// This is the abstract form of a single particle operator
// implemented in functional paradigma
template<typename T>
using SPOF = function<void(const PrimitiveBasis&, Tensor<T>&, const Tensor<T>&)>;

template<typename T>
class SingleParticleOperatorFunction
	: public SPO<T>
	/**
	 * \class singleparticleoperatorfunction
	 * \ingroup operators
	 * \brief this class allows to use functions as spos.
	 */
{
public:
	SingleParticleOperatorFunction() = default;

	SingleParticleOperatorFunction(const SPOF<T> h)
		: h_(move(h)) {}

	~SingleParticleOperatorFunction() = default;

	void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff) const override;

protected:
	SPOF<T> h_;
};

template<typename T>
using SPOf = SingleParticleOperatorFunction<T>;

typedef SingleParticleOperatorFunction<complex<double>> SPOfcd;

