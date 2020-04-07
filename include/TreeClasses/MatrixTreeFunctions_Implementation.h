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
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket, const Node& node) {
		if (!node.isBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				const Node& child = node.child(k);
				Ket = MatrixTensor(S[child], Ket, k);
			}
		}
		size_t last = node.shape().lastIdx();
		Contraction(S[node], Bra, Ket, last);

	}

	template<typename T>
	void DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree) {
		for (const Node& node : tree) {
			DotProductLocal(S, Psi[node], Chi[node], node);
		}
	}

	template<typename T>
	MatrixTree<T> DotProduct(const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		MatrixTree<T> S(tree);
		DotProduct(S, Bra, Ket, tree);
		return S;
	}

////////////////////////////////////////////////////////////////////////
/// General Contraction for Tensor Trees
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void ContractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
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
					Ket = MatrixTensor(S[child], Ket, k);
				}
			}
		}

		Ket = multStateAB(Rho[parent], Ket);

		Contraction(Rho[node], Bra, Ket, child_idx);
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const Tree& tree, const MatrixTree<T> *S_opt = nullptr) {
		if (S_opt != nullptr) {
			assert(S_opt->size() == tree.nNodes());

		}
		assert(Rho.size() == tree.nNodes());
		assert(Psi.size() == tree.nNodes());
		assert(Chi.size() == tree.nNodes());

		// @TODO: Use reverse iterator
		for (auto it = tree.rbegin(); it != tree.rend(); it++) {
			const Node& node = *it;
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				ContractionLocal(Rho, Psi[parent], Chi[parent], node, S_opt);
			} else {
				Rho[node] = IdentityMatrix<T>(node.shape().lastDimension());
			}
		}
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree) {
		Contraction(Rho, Psi, Chi, tree, &S);
	}

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree) {
		MatrixTree<T> Rho(tree);
		Contraction(Rho, Psi, Chi, S, tree);
		return Rho;
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const Tree& tree, bool orthogonal) {
		if (orthogonal) {
			Contraction(Rho, Psi, Psi, tree);
		} else {
			MatrixTree<T> S = DotProduct(Psi, Psi, tree);
			Contraction(Rho, Psi, Psi, S, tree);
		}
	}

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const Tree& tree, bool orthogonal) {
		MatrixTree<T> Rho(tree);
		Contraction(Rho, Psi, tree, orthogonal);
		return Rho;
	}

}

#endif //MATRIXTREE_IMPLEMENTATION_H
