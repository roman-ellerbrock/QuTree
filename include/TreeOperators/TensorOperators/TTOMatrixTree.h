//
// Created by Roman Ellerbrock on 8/4/21.
//

#ifndef TTOMATRIXTREE_H
#define TTOMATRIXTREE_H
#include "TensorTreeOperator.h"
#include "TreeClasses/MatrixTree.h"
#include "TreeOperators/SumOfProductsOperator.h"

template <typename T>
void toTensor(Tensor<T>& A, const MLO<T>& M, size_t part, const Leaf& leaf);

template <typename T>
Tensor<T> toTensor(const SOP<T>& S, const Leaf& leaf);

template <typename T>
T prodMk(const vector<size_t>& idx, const vector<Matrix<T>>& Mk, size_t l, int skip = -1);

template <typename T>
class TTOMatrixTree: public MatrixTree<T> {
	/**
	 * \brief this class calculates the matrix representation required to contract SOPs into a TTNO
	 * \ingroup TTNO
	 */
public:
	using MatrixTreed::NodeAttribute<Matrix<T>>::attributes_;
	using MatrixTreed::NodeAttribute<Matrix<T>>::operator[];

	TTOMatrixTree(const SOP<T>& S, const Tree& tree) {
		size_t npart = S.size();
		attributes_.clear();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrix<T>(ntensor, npart));
		}
	}

	~TTOMatrixTree() = default;

	vector<Matrix<T>> gatherMk(const Node& node) const {
		vector<Matrix<T>> Mk;
		if (node.isBottomlayer()) { cerr << "Bottomlayer does not have active matrices.\n"; getchar(); }
		if (node.isBottomlayer()) { return Mk; }
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			Mk.emplace_back((*this)[child]);
		}
		return Mk;
	}

	void representLayer(const TensorTreeOperator<T>& A, const SOP<T>& S, const Node& node) {
		auto& mrep = (*this)[node];
		mrep.zero();
		if (node.isBottomlayer()) {
			Tensor<T> C = toTensor(S, node.getLeaf());
			mrep = A[node].dotProduct(C);
		} else {
			/// collect all underlying matrices
			/// contribution from each matrix
			const Tensor<T>& B = A[node];
			const TensorShape& shape = B.shape();
			const auto& Mk = gatherMk(node);
			for (size_t l = 0; l < S.size(); ++l) {
				for (size_t I = 0; I < shape.totalDimension(); ++I) {
					auto idx = indexMapping(I, shape);
					size_t i0 = idx[shape.lastIdx()];
					T factor = prodMk(idx, Mk, l);
					mrep(i0, l) += B(I) * factor;
				}
			}
		}
	}

	void represent(const TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& tree) {
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

typedef TTOMatrixTree<double> TTOMatrixTreed;
typedef TTOMatrixTree<complex<double>> TTOMatrixTreecd;

#endif //TTOMATRIXTREE_H
