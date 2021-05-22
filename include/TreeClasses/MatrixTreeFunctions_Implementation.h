//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef MATRIXTREE_IMPLEMENTATION_H
#define MATRIXTREE_IMPLEMENTATION_H
#include "MatrixTreeFunctions.h"

namespace TreeFunctions {
////////////////////////////////////////////////////////////////////////
/// General DotProduct for Tensor Trees
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void dotProductLocal(MatrixTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket, const Node& node) {
		if (!node.isBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				const Node& child = node.child(k);
				Ket = matrixTensor(S[child], Ket, k);
			}
		}
		size_t last = node.shape().lastIdx();
		contraction(S[node], Bra, Ket, last);
		if (S[node].dim1() != S[node].dim2()) {
			cout << "Ket:\n";
			Ket.print();
			cout << "S:\n";
			S[node].print();
			@TODO: Something is wrong in the contraction routine for asymmetric inputs
		}

	}

	template<typename T>
	void dotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree) {
		for (const Node& node : tree) {
			dotProductLocal(S, Psi[node], Chi[node], node);
		}
	}

	template<typename T>
	MatrixTree<T> dotProduct(const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		MatrixTree<T> S(tree);
		dotProduct(S, Bra, Ket, tree);
		return S;
	}

////////////////////////////////////////////////////////////////////////
/// General Contraction for Tensor Trees
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void contractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
		const Node& node, const MatrixTree<T> *S_opt) {
		assert(!node.isToplayer());

		const Node& parent = node.parent();
		auto child_idx = (size_t) node.childIdx();

		/// Optional Overlap matrix
		if (S_opt != nullptr) {
			const MatrixTree<T>& S = *S_opt;
			for (size_t k = 0; k < parent.nChildren(); ++k) {
				if (k != child_idx) {
					const Node& child = parent.child(k);
					Ket = matrixTensor(S[child], Ket, k);
				}
			}
		}

		Ket = multStateAB(Rho[parent], Ket);

		contraction(Rho[node], Bra, Ket, child_idx);
	}

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const Tree& tree, const MatrixTree<T> *S_opt = nullptr) {
		if (S_opt != nullptr) {
			assert(S_opt->size() == tree.nNodes());
		}
		assert(Rho.size() == tree.nNodes());
		assert(Psi.size() == tree.nNodes());
		assert(Chi.size() == tree.nNodes());

		for (auto it = tree.rbegin(); it != tree.rend(); it++) {
			const Node& node = *it;
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				contractionLocal(Rho, Psi[parent], Chi[parent], node, S_opt);
			} else {
				Rho[node] = identityMatrix<T>(node.shape().lastDimension());
			}
		}
	}

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree) {
		contraction(Rho, Psi, Chi, tree, &S);
	}

	template<typename T>
	MatrixTree<T> contraction(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree) {
		MatrixTree<T> Rho(tree);
		contraction(Rho, Psi, Chi, S, tree);
		return Rho;
	}

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const Tree& tree, bool orthogonal) {
		if (orthogonal) {
			contraction(Rho, Psi, Psi, tree);
		} else {
			MatrixTree<T> S = dotProduct(Psi, Psi, tree);
			contraction(Rho, Psi, Psi, S, tree);
		}
	}

	template<typename T>
	MatrixTree<T> contraction(const TensorTree<T>& Psi, const Tree& tree, bool orthogonal) {
		MatrixTree<T> Rho(tree);
		contraction(Rho, Psi, tree, orthogonal);
		return Rho;
	}

}

#endif //MATRIXTREE_IMPLEMENTATION_H
