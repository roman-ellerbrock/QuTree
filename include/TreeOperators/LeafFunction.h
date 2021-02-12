#pragma once
#include "LeafOperator.h"

template<typename T>
using LeafFun = function<void(const LeafInterface&, Tensor<T>&, const Tensor<T>&)>;

template<typename T>
using LeafFunPair = pair<LeafFun<T>, LeafFun<T>>;

typedef LeafFunPair<complex<double>> LeafFunPaircd;
typedef LeafFunPair<double> LeafFunPaird;

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

	explicit LeafFunction(const LeafFun<T> h)
		: h_(move(h)) {}

	~LeafFunction() = default;

	void apply(const LeafInterface& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff) const override;

protected:
	LeafFun<T> h_;
};

typedef LeafFun<complex<double>> LeafFuncd;
typedef LeafFunction<complex<double>> LeafFunctioncd;

