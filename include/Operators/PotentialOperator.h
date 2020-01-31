#pragma once
#include "Core/stdafx.h"

/**
 * \class PotentialOperator
 * \brief This Operator tells if a potential operator is present.
 *
 * This operator tells that a potential acts is present. It also manages
 * the number of degrees of freedom used in the PES and the state used
 * in the PES.
 * */


class PotentialOperator
	/**
	 * \class PotentialOperator
	 * \ingroup Operators
	 * \brief This Operator marks a potential operator in a MultiParticleOperator
	 *
	 * The class is only relevant when working with quadrature methods on tensor trees.
	 */
{
public:
	PotentialOperator()
	  :f(0),state(0){}

	PotentialOperator(size_t f_, size_t state_)
	:f(f_), state(state_){}

	~PotentialOperator() = default;

	size_t F()const {
		return f;
	}

	size_t State()const {
		return state;
	}
	
protected:
	size_t f;
	size_t state;

};


