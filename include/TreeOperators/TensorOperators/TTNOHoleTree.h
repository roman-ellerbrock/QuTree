//
// Created by Roman Ellerbrock on 8/6/21.
//

#ifndef TTNOHOLETREE_H
#define TTNOHOLETREE_H
#include "TTNOMatrixTree.h"

class TTNOHoleTree: public MatrixTreed {
	using MatrixTreed::NodeAttribute<Matrix<double>>::
	attributes_;
public:
	TTNOHoleTree() = default;
	~TTNOHoleTree() = default;

	TTNOHoleTree(const SOPd& S, const Tree& tree) {
		size_t npart = S.size();
		attributes_.clear();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrixd(ntensor, npart));
		}
	}

	void gatherMk(vector<Matrixd>& Mk, const Node& node) const {
		if (!node.isToplayer()) {
			Mk.push_back((*this)[node.parent()]);
		}
	}

	void represent(const TensorOperatorTree& A, const TTNOMatrixTree& M, const Tree& tree) {
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

			for (size_t l = 0; l < mrep.dim2(); ++l) {
				for (size_t I = 0; I < shape.totalDimension(); ++I) {
					auto is = indexMapping(I, shape);
					size_t ik = is[skip];
					double factor = prodMk(is, ms, l, (int) skip);
					mrep(ik, l) += B(I) * factor;
				}
			}
		}
	}
};


#endif //TTNOHOLETREE_H
