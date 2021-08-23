//
// Created by Roman Ellerbrock on 8/6/21.
//

#ifndef TTOHOLETREE_H
#define TTOHOLETREE_H
#include "TTOMatrixTree.h"

class TTOHoleTree: public MatrixTreed {
	/**
	 * \brief This class is the contraction required to contract a SOP into a TTNO.
	 * \ingroup TTNO
	 */
	using MatrixTreed::NodeAttribute<Matrix<double>>::attributes_;
public:
	TTOHoleTree() = default;
	~TTOHoleTree() = default;

	TTOHoleTree(const SOPd& S, const Tree& tree) {
		size_t npart = S.size();
		attributes_.clear();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrixd(ntensor, npart));
		}
	}

	/**
	 *
	 * @param Mk
	 * @param parent of node at which calculation is performed (NOT node itself!)
	 */
	void gatherMk(vector<Matrixd>& Mk, const Node& parent) const {
		if (!parent.isToplayer()) {
			Mk.push_back((*this)[parent]);
		}
	}

	void represent(const TensorTreeOperator& A, const TTOMatrixTree& M, const Tree& tree) {
		for (int no = tree.nNodes() - 2; no >= 0; no--) {
			const Node& node = tree.getNode(no);
			const Node& parent = node.parent();
			size_t skip = node.childIdx();
			/// collect all underlying matrices
			vector<Matrixd> ms = M.gatherMk(parent);
			gatherMk(ms, parent);

			/// contribution from each matrix
			const Tensord& B = A[parent];
			const TensorShape& shape = B.shape();
			auto& mrep = (*this)[node];
			mrep.zero();

			for (size_t l = 0; l < mrep.dim2(); ++l) {
				for (size_t I = 0; I < shape.totalDimension(); ++I) {
					auto idx = indexMapping(I, shape);
					size_t ik = idx[skip];
					double factor = prodMk(idx, ms, l, (int) skip);
					mrep(ik, l) += B(I) * factor;
				}
			}
		}
	}
};


#endif //TTOHOLETREE_H