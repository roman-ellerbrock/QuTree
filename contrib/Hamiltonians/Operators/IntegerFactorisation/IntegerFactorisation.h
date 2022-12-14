//
// Created by Roman Ellerbrock on 2019-07-12.
//

#ifndef INTEGERFACTORISATION_H
#define INTEGERFACTORISATION_H
#include "SumOfProductsOperator.h"
#include "Pauli.h"

class IntegerFactorisation: public SOP {
	/*!
	 * \brief This is an implementation of the number partitioning problem as an Ising model
	 *
	 * The number partitioning problem is the problem of partitioning a set of numbers into two subsets
	 * with equal sums, e.g. A = (1, 2, 3, 4) -> G1 = (1, 4), G2 = (2, 3) (both elements in G1 and G2
	 * add up to 5.
	 * The problem can be mapped to the hamiltonian H = sum_i (n(i) * s_z(i))^2, where n(i) is the
	 * i-th number in the set A and s_z(i) is the pauli z-matrix acting on qubit i (see Ref. [1]).
	 *
	 *
	 * [1] DOI: 10.3389/fphy.2014.00005
	 */
public:
	/// Constructor
	explicit IntegerFactorisation(const mctdhBasis& basis) { Initialize(basis); }
	IntegerFactorisation(const mctdhBasis& basis, size_t F);


	/// Destructor
	~IntegerFactorisation() = default;

	/// This routine initializes the SOP operator which is called implicitely by Initialize routine
	void SpecialInitialize(const mctdhBasis& basis, const size_t F);
};

#endif //NUMBERPARTITIONING_H
