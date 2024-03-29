//
// Created by Roman Ellerbrock on 8/6/21.
//

#ifndef TTOHOLETREE_H
#define TTOHOLETREE_H
#include "TTOMatrixTree.h"

template<typename T>
class TTOHoleTree: public NodeAttribute<Matrix<T>> {
	/**
	 * \brief This class is the contraction required to contract a SOP into a TTNO.
	 * \ingroup TTNO
	 */
public:
	using NodeAttribute<Matrix<T>>::attributes_;
	TTOHoleTree() = default;
	~TTOHoleTree() = default;

	TTOHoleTree(const SOP<T>& S, const Tree& tree) {
		size_t npart = S.size();
		attributes_.clear();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrix<T>(ntensor, npart));
		}
	}

	/**
	 *
	 * @param Mk
	 * @param parent of node at which calculation is performed (NOT node itself!)
	 */
	void gatherMk(vector<Matrix<T>>& Mk, const Node& parent) const {
		if (!parent.isToplayer()) {
			Mk.push_back((*this)[parent]);
		}
	}

	void represent(const TensorTreeOperator<T>& A, const TTOMatrixTree<T>& M, const Tree& tree) {
		for (int no = tree.nNodes() - 2; no >= 0; no--) {
			const Node& node = tree.getNode(no);
			const Node& parent = node.parent();
			size_t skip = node.childIdx();
			/// collect all underlying matrices
			vector<Matrix<T>> ms = M.gatherMk(parent);
			gatherMk(ms, parent);

			/// contribution from each matrix
			const Tensor<T>& B = A[parent];
			const TensorShape& shape = B.shape();
			auto& mrep = (*this)[node];
			mrep.zero();

			auto idx = indexMapping(0, shape);
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				indexMapping(idx, I, shape);
				size_t ik = idx[skip];
				if (ik != 0) { continue; }
				for (size_t l = 0; l < mrep.dim2(); ++l) {
					T factor = prodMk(idx, ms, l, (int) skip);
					for (size_t ik = 0; ik < shape[skip]; ++ik) {
						size_t J = I + shape.before(skip) * ik;
						mrep(ik, l) += B(J) * factor;
					}
				}
			}
		}
	}
};

typedef TTOHoleTree<double> TTOHoleTreed;

typedef TTOHoleTree<complex<double>> TTOHoleTreecd;

#endif //TTOHOLETREE_H
