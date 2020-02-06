#pragma once
#include "LeafOperator.h"

// This is the abstract form of a single particle operator
// implemented in functional paradigma
template<typename T>
using SPOF = function<void(const PrimitiveBasis&, Tensor<T>&, const Tensor<T>&)>;

template<typename T>
class LeafFunction
	: public LeafOperator<T>
	/**
	 * \class LeafFunction
	 * \ingroup Operators
	 * \brief This class allows to use functions as LeafOperators.
	 */
{
public:
	LeafFunction() = default;

	explicit LeafFunction(const SPOF<T> h)
		: h_(move(h)) {}

	~LeafFunction() = default;

	void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff) const override;

protected:
	SPOF<T> h_;
};

template<typename T>
using SPOf = LeafFunction<T>;

typedef LeafFunction<complex<double>> SPOfcd;

