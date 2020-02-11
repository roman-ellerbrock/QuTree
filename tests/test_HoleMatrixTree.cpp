//
// Created by Roman Ellerbrock on 2/10/20.
//
#include "UnitTest++/UnitTest++.h"
#include "IOTree.h"

SUITE(HoleMatrixTree) {
	TEST(IO) {
		mt19937 gen(1993);
		TTBasis basis(14, 4, 2);
		TensorTreecd Psi(basis, gen, false);
		IOTree::Occupancy(Psi, basis);
		HoleMatrixTreecd Rho(Psi, basis);
		IOTree::Leafs(Psi, basis);
	}


}


