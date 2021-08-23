//
// Created by Roman Ellerbrock on 8/15/21.
//

#ifndef TTOCONTRACTION_H
#define TTOCONTRACTION_H
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeOperators/TensorOperators/TTOrepresentation.h"

class TTOcontraction: public NodeAttribute<vector<Matrixd>> {
	/**
	 * \brief This class calculates the contraction of a TTNO representation
	 * \ingroup TTNO
	 *
	 */
public:
	TTOcontraction(const Tree& tree, const Tree& optree) {
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
				const Node& opnode = optree.getNode(node.address());
				const
				size_t nSPF = node.shape().lastDimension();
				size_t nSPO = opnode.shape().lastDimension();
				Matrixd h(nSPF, nSPF);
				vector<Matrixd> hs(nSPO, h);
				attributes_.emplace_back(hs);
			}
		}
	}

	~TTOcontraction() = default;

	Tensord applyMatrices(Tensord A, const vector<size_t>& ls, const TTOrepresentation& Hrep,
		const Node& parent, size_t hole) {
		for (size_t k = 0; k < parent.nChildren(); ++k) {
			if (k == hole) { continue; }
			const Node& child = parent.child(k);
			const vector<Matrixd>& hrep = Hrep[child];
			A = matrixTensor(hrep[ls[k]], A, k);
		}
		if (!parent.isToplayer()) {
			const vector<Matrixd> hholes = (*this)[parent];
			size_t l0 = ls[parent.parentIdx()];
			const Matrixd& hhole = hholes[l0];
			A = matrixTensor(hhole, A, parent.parentIdx());
		}
		return A;
	}

	void calculateLayer(const Tensord& Psi, const TensorTreeOperator& H,
		const TTOrepresentation& Hrep, const Tensord& Chi, const Node& opnode) {

		if (opnode.isToplayer()) {
			cerr << "Error: child may not be root.\n";
			exit(1);
		}
		const Node& parent = opnode.parent();
		const Tensord& B = H[parent];
		const TensorShape& opshape = B.shape();
		size_t holeidx = opnode.childIdx();

		vector<Matrixd>& hs = (*this)[opnode];
		for (Matrixd& h : hs) {
			h.zero();
		}
		for (size_t L = 0; L < opshape.totalDimension(); ++L) {
			auto ls = indexMapping(L, opshape);
			auto hChi = applyMatrices(Chi, ls, Hrep, parent, holeidx);
			auto hk = contraction(Psi, hChi, holeidx);
			hs[ls[holeidx]] += B[L] * hk;
		}
	}

	void calculate(const TensorTreed& Psi, const TensorTreeOperator& H,
		const TTOrepresentation& Hrep, const TensorTreed& Chi, const Tree& optree) {
		for (int l = optree.nNodes() - 2; l >= 0; l--) {
			const Node& node = optree.getNode(l);
			if (node.isToplayer()) { continue; }
			const Node& parent = node.parent();
			calculateLayer(Psi[parent], H, Hrep, Chi[parent], node);
		}
	}

	void print(const Tree& tree) {
		for (const Node& node : tree) {
			if (node.isToplayer()) { continue; }
			const auto& hs = (*this)[node];
			node.info();
			for (size_t l = 0; l < hs.size(); ++l) {
				cout << "l = " << l << endl;
				hs[l].print();
			}
		}
	}
};


#endif //TTOCONTRACTION_H

