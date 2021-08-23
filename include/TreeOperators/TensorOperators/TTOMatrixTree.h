//
// Created by Roman Ellerbrock on 8/4/21.
//

#ifndef TTOMATRIXTREE_H
#define TTOMATRIXTREE_H
#include "TensorTreeOperator.h"
#include "TreeClasses/MatrixTree.h"
#include "TreeOperators/SumOfProductsOperator.h"

void toTensor(Tensord& A, const MLOd& M, size_t part, const Leaf& leaf);
Tensord toTensor(const SOPd& S, const Leaf& leaf);

double prodMk(const vector<size_t>& idx, const vector<Matrixd>& Mk, size_t l, int skip = -1);

class TTOMatrixTree: public MatrixTreed {
	/**
	 * \brief this class calculates the matrix representation required to contract SOPs into a TTNO
	 * \ingroup TTNO
	 */
	using MatrixTreed::NodeAttribute<Matrix<double>>::attributes_;
public:
	using MatrixTreed::NodeAttribute<Matrix<double>>::operator[];

	TTOMatrixTree(const SOPd& S, const Tree& tree) {
		size_t npart = S.size();
		attributes_.clear();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrixd(ntensor, npart));
		}
	}

	~TTOMatrixTree() = default;

	vector<Matrixd> gatherMk(const Node& node) const {
		vector<Matrixd> Mk;
		if (node.isBottomlayer()) { cerr << "Bottomlayer does not have active matrices.\n"; getchar(); }
		if (node.isBottomlayer()) { return Mk; }
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			Mk.emplace_back((*this)[child]);
		}
		return Mk;
	}

	void representLayer(const TensorTreeOperator& A, const SOPd& S, const Node& node) {
		auto& mrep = (*this)[node];
		mrep.zero();
		if (node.isBottomlayer()) {
			Tensord C = toTensor(S, node.getLeaf());
			mrep = A[node].dotProduct(C);
		} else {
			/// collect all underlying matrices
			/// contribution from each matrix
			const Tensord& B = A[node];
			const TensorShape& shape = B.shape();
			const auto& Mk = gatherMk(node);
			for (size_t l = 0; l < S.size(); ++l) {
				for (size_t I = 0; I < shape.totalDimension(); ++I) {
					auto idx = indexMapping(I, shape);
					size_t i0 = idx[shape.lastIdx()];
					double factor = prodMk(idx, Mk, l);
					mrep(i0, l) += B(I) * factor;
				}
			}
		}
	}

	void represent(const TensorTreeOperator& A, const SOPd& S, const Tree& tree) {
		for (const Node& node : tree) {
			representLayer(A, S, node);
		}
	}

	void print(const Tree& tree) {
		for (const Node& node : tree) {
			node.info();
			(*this)[node].print();
		}
	}

};


#endif //TTOMATRIXTREE_H