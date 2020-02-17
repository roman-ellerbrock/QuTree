//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef MATRIXTREE_IMPLEMENTATION_H
#define MATRIXTREE_IMPLEMENTATION_H
#include "MatrixTreeFunctions.h"

namespace MatrixTreeFunctions {
////////////////////////////////////////////////////////////////////////
/// General DotProduct for Tensor Trees
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket, const Node& node) {
		if (!node.IsBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				const Node& child = node.Down(k);
				Ket = multAB(S[child], Ket, k);
			}
		}
		size_t last = node.TDim().GetLastIdx();
		mHoleProduct(S[node], Bra, Ket, last);
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
		assert(!node.IsToplayer());

		const Node& parent = node.Up();
		auto child_idx = (size_t) node.ChildIdx();

		/// Optional Overlap matrix
		if (S_opt != nullptr) {
			const MatrixTree<T>& S = *S_opt;
			for (size_t k = 0; k < parent.nChildren(); ++k) {
				if (k != child_idx) {
					const Node& child = parent.Down(k);
					Ket = multAB(S[child], Ket, k);
				}
			}
		}

		if (!parent.IsToplayer()) {
			Ket = multStateAB(Rho[parent], Ket);
		}

		mHoleProduct(Rho[node], Bra, Ket, child_idx);
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
		for (int n = (int) tree.nNodes() - 2; n >= 0; --n) {
			const Node& node = tree.GetNode(n);
			const Node& parent = node.Up();
			ContractionLocal(Rho, Psi[parent], Chi[parent], node, S_opt);
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