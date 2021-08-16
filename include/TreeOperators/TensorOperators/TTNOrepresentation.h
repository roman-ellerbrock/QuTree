//
// Created by Roman Ellerbrock on 8/15/21.
//

#ifndef TTNOREPRESENTATION_H
#define TTNOREPRESENTATION_H
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/TensorOperators/TensorOperatorTree.h"

class TTNOrepresentation : public NodeAttribute<vector<Matrixd>> {
public:
	TTNOrepresentation(const Tree& tree, const Tree& optree) {
		attributes_.clear();
		for (const Node& node : tree) {
			const Node& opnode = optree.getNode(node.address());
			size_t nSPF = node.shape().lastDimension();
			size_t nSPO = opnode.shape().lastDimension();
			TensorShape shape({nSPF, nSPF, nSPO});
			Matrixd h(nSPF, nSPF);
			vector<Matrixd> hs(nSPO, h);
			attributes_.emplace_back(hs);
		}
	}

	~TTNOrepresentation() = default;

	Tensord applyMatrices(Tensord A, const Tensord& B, const size_t l,
		const Leaf& leaf) {
		size_t dim = leaf.dim();
		Matrixd h(dim, dim);
		for (size_t I = 0; I < dim*dim; ++I) {
			h[I] = B(I, l);
		}
		A = matrixTensor(h, A, 0);
		return A;
	}

	Tensord applyMatrices(Tensord A, const vector<size_t>& ls,
		const Node& opnode) {
		const TensorShape& opshape = opnode.shape();
		for (size_t k = 0; k < opnode.nChildren(); ++k) {
			const Node& child = opnode.child(k);
			const vector<Matrixd>& hs = (*this)[child];
			A = matrixTensor(hs[ls[k]], A, k);
		}
		return A;
	}

	void calculateLayer(const Tensord& Psi, const TensorOperatorTree& H,
		const Tensord& Chi, const Node& opnode) {

		const TensorShape& shape = Psi.shape();
		const Tensord& B = H[opnode];
		const TensorShape& opshape = B.shape();

		vector<Matrixd>& hs = (*this)[opnode];
		if (opnode.isBottomlayer()) {
			const Leaf& leaf = opnode.getLeaf();
			leaf.info();
			for (size_t l = 0; l < opshape.lastDimension(); ++l) {
				auto hChi = applyMatrices(Chi, H[opnode], l, leaf);
				auto hij = Psi.dotProduct(hChi);
				hs[l] = hij;
			}
		} else {
			for (size_t L = 0; L < opshape.totalDimension(); ++L) {
				auto ls = indexMapping(L, opshape);
				auto hChi = applyMatrices(Chi, ls, opnode);
				auto hij = Psi.dotProduct(hChi);
				hs[ls[opnode.parentIdx()]] += B[L] * hij;
			}
		}
	}

	void calculate(const TensorTreed& Psi, const TensorOperatorTree& H, const TensorTreed& Chi,
		const Tree& optree) {
		for (const Node& node : optree) {
			calculateLayer(Psi[node], H, Chi[node], node);
		}
	}

	void print(const Tree& tree) {
		for (const Node& node : tree) {
			const auto& hs = (*this)[node];
			node.info();
			for (size_t l = 0; l < hs.size(); ++l) {
				cout << "l = " << l << endl;
				hs[l].print();
			}
		}
	}

};


#endif //TTNOREPRESENTATION_H
