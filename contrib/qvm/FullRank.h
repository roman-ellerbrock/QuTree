//
// Created by Roman Ellerbrock on 4/14/20.
//

#ifndef FULLRANK_H
#define FULLRANK_H
#include "TreeShape/Tree.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SOPVector.h"
#include "Measurements.h"

namespace FullRank {
	typedef Tensorcd Wavefunction;

    Matrixcd toMatrix(const LeafOperatorcd& h, size_t mode, const Tree& tree);
	Wavefunction applyOperator(Wavefunction Psi, const SOPVectorcd& S, const Tree& tree);

	Wavefunction initialize(const Tree& tree);

	void print(const Wavefunction& Psi);

	complex<double> probability(const vector<size_t>& configuration, const Wavefunction& Psi);
	double probability(const FullRank::Wavefunction& Psi, size_t idx);
	void normalize(Wavefunction& Psi);
}

#endif //FULLRANK_H
