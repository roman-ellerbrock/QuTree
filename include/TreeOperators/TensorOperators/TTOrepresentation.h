//
// Created by Roman Ellerbrock on 8/15/21.
//

#ifndef TTOREPRESENTATION_H
#define TTOREPRESENTATION_H
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"

template <typename T>
class TTOrepresentation : public NodeAttribute<vector<Matrix<T>>> {
	/**
	 * \brief Representation of TTNO in TTN basis
	 * \ingroup TTNO
	 *
	 * This class calculates the matrix representations of TTNO
	 * that are required for applying TTNOs to wavefunctions.
	 */
	 using NodeAttribute<vector<Matrix<T>>>::attributes_;
public:

	TTOrepresentation(const Tree& tree, const Tree& optree) {
		attributes_.clear();
		for (const Node& node : tree) {
			const Node& opnode = optree.getNode(node.address());
			size_t nSPF = node.shape().lastDimension();
			size_t nSPO = opnode.shape().lastDimension();
			Matrix<T> h(nSPF, nSPF);
			vector<Matrix<T>> hs(nSPO, h);
			attributes_.emplace_back(hs);
		}
	}

	~TTOrepresentation() = default;

	Tensor<T> applyMatrices(Tensor<T> A, const Tensor<T>& B, const size_t l,
		const Leaf& leaf) {
		size_t dim = leaf.dim();
		Matrix<T> h(dim, dim);
		for (size_t I = 0; I < dim*dim; ++I) {
			h[I] = B(I, l);
		}
		A = matrixTensor(h, A, 0);
		return A;
	}

	Tensor<T> applyMatrices(Tensor<T> A, const vector<size_t>& ls,
		const Node& opnode) {
		for (size_t k = 0; k < opnode.nChildren(); ++k) {
			const Node& child = opnode.child(k);
			const vector<Matrix<T>>& hs = (*this)[child];
			A = matrixTensor(hs[ls[k]], A, k);
		}
		return A;
	}

	void calculateLayer(const Tensor<T>& Psi, const TensorTreeOperator<T>& H,
		const Tensor<T>& Chi, const Node& opnode) {

		const TensorShape& shape = Psi.shape();
		const Tensor<T>& B = H[opnode];
		const TensorShape& opshape = B.shape();

		vector<Matrix<T>>& hs = (*this)[opnode];
		for (Matrix<T>& h : hs) {
			h.zero();
		}
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

	void calculate(const TensorTree<T>& Psi, const TensorTreeOperator<T>& H, const TensorTree<T>& Chi,
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

typedef TTOrepresentation<complex<double>> TTOrepresentationcd;
typedef TTOrepresentation<double> TTOrepresentationd;

#endif //TTOREPRESENTATION_H
