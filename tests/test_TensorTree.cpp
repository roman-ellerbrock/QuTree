//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "Tree/TensorTree.h"
#include "Tree/TensorTree_Implementation.h"
#include "TensorTreeBasis/TTBasisFactory.h"

SUITE (TensorTree) {

	TEST (TensorTree_FILE_IO) {
		TensorTreeBasis basis(12, 2, 4);
		TensorTreecd T(basis);
		T.Write("TT.tmp.dat");
		TensorTreecd Q(basis);
		ifstream is("TT.tmp.dat");
		is >> Q;
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_RandomGenerate) {
		TensorTreeBasis basis(12, 2, 2);
		TensorTreecd T(basis);
		mt19937 gen(2468);
		T.Generate(basis, gen, false);
		string filename("TT.RNG.dat");
		T.Write(filename);
		TensorTreecd Q(filename);
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_Train) {
		size_t nLeaves = 12;
		auto basis = TTBasisFactory::TensorTrain(nLeaves, 4, 2, 6);
			CHECK_EQUAL(2 * nLeaves - 1, basis.nNodes());
		mt19937 gen(2468);
		TensorTreecd Psi(basis, gen);
			CHECK_EQUAL(basis.nNodes(), Psi.size());
	}

}

