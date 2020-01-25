//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TensorTreeBasis.h"
#include "TensorTree.h"
#include "TensorTree_Implementation.h"
#include "DenseOverlap.h"
#include "HoleOverlap.h"

SUITE (Tree) {

	TEST (TensorTreeBasis_Generator) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);

		for (size_t n_modes = 2; n_modes < 18; ++n_modes) {
			TensorTreeBasis basis(n_modes, n_leaf, n_node);

			for (const Node& node : basis) {
				const TensorDim& tdim = node.TDim();
				if (node.IsToplayer()) {
						CHECK_EQUAL(tdim_top, tdim);
				} else if (node.IsBottomlayer()) {
						CHECK_EQUAL(tdim_bottom, tdim);
				} else {
						CHECK_EQUAL(tdim_upper, tdim);
				}
			}
		}
	}

	TEST (TensorTreeBasis_FileIO) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);
		size_t n_modes = 13;

		TensorTreeBasis basis(n_modes, n_leaf, n_node);
		{
			ofstream os("TTBasis.IO.tmp.dat");
			basis.Write(os);
			os.close();
		}
		TTBasis basis2("TTBasis.IO.tmp.dat");
			CHECK_EQUAL(basis2.nNodes(), basis.nNodes());
	}

	TEST (TensorTreeBasis_Reindexing) {
		size_t n_modes = 9;
		TensorTreeBasis basis(n_modes, 2, 4);

		map<size_t, size_t> Map;
		for (size_t k = 0; k < n_modes; ++k) {
			Map[k] = n_modes - 1 - k;
		}
		basis.ReindexLeafModes(Map);

		size_t k = 0;
		for (const Node& node : basis) {
			if (node.IsBottomlayer()) {
				const Leaf& leaf = node.PhysCoord();
					CHECK_EQUAL(k++, leaf.Mode());
			}
		}
	}

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
		T.Generate(basis, gen,false);
		string filename("TT.RNG.tmp.dat");
		T.Write(filename);
		TensorTreecd Q(filename);
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_Overlap) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.tmp.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		const FactorMatrixcd& s = S.Get();
			CHECK_CLOSE(1., abs(s[0]), 1e-14);
	}

	TEST (DenseOverlap_IO) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.tmp.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		string file2("DO.tmp.dat");
		{
			ofstream os(file2);
			S.Write(os);
			os.close();
		}
		DenseOverlapcd Q(file2);
			CHECK_EQUAL(S.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(S[node], Q[node]);
		}
	}

	TEST (TensorTree_HoleOverlap) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.tmp.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		HoleOverlapcd Rho(T, T, S, basis);
	}

	TEST (HoleOverlap_IO) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.tmp.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		HoleOverlapcd Rho(T, T, S, basis);
		string file2("HO.tmp.dat");
		{
			ofstream os(file2);
			Rho.Write(os);
			os.close();
		}
		HoleOverlapcd Qh(file2);
			CHECK_EQUAL(Rho.size(), Qh.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(Rho[node], Qh[node]);
		}
	}
}

