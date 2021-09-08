//
// Created by Roman Ellerbrock on 8/15/21.
//

#ifndef TTOCONTRACTION_H
#define TTOCONTRACTION_H
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeOperators/TensorOperators/TTOrepresentation.h"

template <typename T>
class TTOcontraction: public NodeAttribute<vector<Matrix<T>>> {
	/**
	 * \brief This class calculates the contraction of a TTNO representation
	 * \ingroup TTNO
	 *
	 */
	using NodeAttribute<vector<Matrix<T>>>::attributes_;
public:
	TTOcontraction() = default;

	TTOcontraction(const Tree& tree, const Tree& optree) {
		initialize(tree, optree);
	}

	void initialize(const Tree& tree, const Tree& optree) {
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
				const Node& opnode = optree.getNode(node.address());
				const
				size_t nSPF = node.shape().lastDimension();
				size_t nSPO = opnode.shape().lastDimension();
				Matrix<T> h(nSPF, nSPF);
				vector<Matrix<T>> hs(nSPO, h);
				attributes_.emplace_back(hs);
			}
		}
	}

	~TTOcontraction() = default;

	[[nodiscard]] Tensor<T> applyMatrices(Tensor<T> A, const vector<size_t>& ls,
		const TensorTreeOperator<T>& H, const TTOrepresentation<T>& Hrep,
		const Node& parent, int hole) const {
/*		for (size_t k = 0; k < parent.nChildren(); ++k) {
			if (k == hole) { continue; }
			const Node& child = parent.child(k);
			const vector<Matrix<T>>& hrep = Hrep[child];
			A = matrixTensor(hrep[ls[k]], A, k);
		}*/
		A = Hrep.applyMatrices(A, H[parent], ls, parent, hole);

		if (!parent.isToplayer() && (hole != parent.parentIdx())) {
			const vector<Matrix<T>>& hholes = (*this)[parent];
			size_t l0 = ls[parent.parentIdx()];
			const Matrix<T>& hhole = hholes[l0];
			A = matrixTensor(hhole, A, parent.parentIdx());
		}
		return A;
	}

	void calculateLayer(const Tensor<T>& Psi, const TensorTreeOperator<T>& H,
		const TTOrepresentation<T>& Hrep, const Tensor<T>& Chi, const Node& opnode) {

		if (opnode.isToplayer()) {
			cerr << "Error: child may not be root.\n";
			exit(1);
		}
		const Node& parent = opnode.parent();
		const Tensor<T>& B = H[parent];
		const TensorShape& opshape = B.shape();
		size_t holeidx = opnode.childIdx();

		vector<Matrix<T>>& hs = (*this)[opnode];
		for (Matrix<T>& h : hs) {
			h.zero();
		}
		auto ls = indexMapping(0, opshape);
		for (size_t L = 0; L < opshape.totalDimension(); ++L) {
			indexMapping(ls, L, opshape);
			auto hChi = applyMatrices(Chi, ls, H, Hrep, parent, holeidx);
			auto hk = contraction(Psi, hChi, holeidx);
//			hs[ls[holeidx]] += B[L] * hk;
			hs[ls[holeidx]] += hk;
		}
	}

	void calculate(const TensorTree<T>& Psi, const TensorTreeOperator<T>& H,
		const TTOrepresentation<T>& Hrep, const TensorTree<T>& Chi, const Tree& optree) {
		for (int l = optree.nNodes() - 2; l >= 0; l--) {
			const Node& node = optree.getNode(l);
			if (node.isToplayer()) { continue; }
			const Node& parent = node.parent();
			calculateLayer(Psi[parent], H, Hrep, Chi[parent], node);
		}
	}

	void print(const Tree& tree) const{
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

template <typename T>
Tensor<T> apply(const Tensor<T>& Phi, const TensorTreeOperator<T>& H, const TTOrepresentation<T>& rep,
	const TTOcontraction<T>& con, const Node& node);

typedef TTOcontraction<double> TTOcontractiond;
typedef TTOcontraction<complex<double>> TTOcontractioncd;

#endif //TTOCONTRACTION_H

