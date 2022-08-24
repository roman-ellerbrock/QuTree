
//
// Created by Roman Ellerbrock on 8/23/22.
//

#include <gtest/gtest.h>
#include "TensorNetwork/ApplySCF.h"
#include "TensorNetwork/contractions.h"
#include "Tree/TreeFactory.h"
#include "Util/GateOperators.h"
#include "TensorNetwork/ExpectationValue.h"
#include "Util/HamiltonianSuite.h"


TEST(ExpectatationValue) {
	Tree tree = balancedTree(4, 10, 4);
	TensorTreecd Psi(tree);
    SOPcd H = CoupledHarmonicOsscillator(4);
    auto exp = expectationValue(Psi, H);
    exp.print();
    getchar();
}