#pragma once
#include "Core/Tensor.h"
#include "PrimitiveBasis.h"

/**
 * \defgroup Operators
 * \brief This group is a bundle of classes that handle
 * operators. System specific operators are NOT included here.
 *
 * This group contains general operator classes that are required to
 * evaluate working equations for tensor tree approaches.
 */

template <typename T>
class LeafOperator
	/**
	 * \class LeafOperator
	 * \ingroup Operators
	 * \brief This class represents a single particle operator acting on a single leaf.
	 */
{
public:
	LeafOperator() = default;
	~LeafOperator() = default;

	virtual void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff)const = 0;
};

typedef LeafOperator<complex<double>> LeafOperatorcd;
