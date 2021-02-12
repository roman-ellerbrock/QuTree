#pragma once
#include "Core/stdafx.h"
#include "TreeOperators/Potential.h"
#include "TreeOperators/CoordinateTransformation.h"

/**
 * \class PotentialOperator
 * \brief This Operator tells if a potential operator is present.
 *
 * This operator tells that a potential acts is present. It also manages
 * the number of degrees of freedom used in the PES and the state_ used
 * in the PES.
 * */


class PotentialOperator
	/**
	 * \class PotentialOperator
	 * \ingroup Operators
	 * \brief This Operator marks a potential operator in a MultiLeafOperator
	 *
	 * The class is only relevant when working with quadrature methods on tensor trees.
	 */
{
public:
	PotentialOperator()
	  : f_(0), state_(0), V_(nullptr), Q_(new CoordinateTransformation()) {}

	PotentialOperator(shared_ptr<Potential> V, size_t f_, size_t state_)
		: f_(f_), state_(state_), V_(V), Q_(new CoordinateTransformation()) {}

	~PotentialOperator() = default;

	size_t F()const {
		return f_;
	}

	size_t State()const {
		return state_;
	}

	shared_ptr<Potential>& V() { return V_; }

	double Evaluate(const Vectord& Xv, size_t part) const {
		Vectord q = Q_->transform(Xv);
		assert(V_);
		return V_->evaluate(q, part);
	}

	shared_ptr<CoordinateTransformation> Q_;

protected:
	size_t f_;
	size_t state_;
	shared_ptr<Potential> V_;
};


