#pragma once
#include "Core/stdafx.h"

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
	  : f_(0), state_(0){}

	PotentialOperator(size_t f_, size_t state_)
	: f_(f_), state_(state_){}

	~PotentialOperator() = default;

	size_t F()const {
		return f_;
	}

	size_t State()const {
		return state_;
	}
	
protected:
	size_t f_;
	size_t state_;

};


