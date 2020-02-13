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
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Phi, Tensor<T> AChi, const Node& node) {
		if (!node.IsBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				const Node& child = node.Down(k);
				const Matrix<T>& spo = S[child];
				AChi = multAB(spo, AChi, k);
			}
		}
		size_t order = node.TDim().GetOrder();
		mHoleProduct(S[node], Phi, AChi, order);
	}

	template<typename T>
	void DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const TTBasis& basis) {
		for (const Node& node : basis) {
			DotProductLocal(S, Psi[node], Chi[node], node);
		}
	}

	template<typename T>
	MatrixTree<T> DotProduct(const TensorTree<T>& Psi, const TensorTree<T>& Chi, const TTBasis& basis) {
		MatrixTree<T> S(basis);
		DotProduct(S, Psi, Chi, basis);
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
		const TTBasis& basis, const MatrixTree<T> *S_opt = nullptr) {
		if (S_opt != nullptr) {
			assert(S_opt->size() == basis.nNodes());
		}
		assert(Rho.size() == basis.nNodes());
		assert(Psi.size() == basis.nNodes());
		assert(Chi.size() == basis.nNodes());

		for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
			const Node& node = basis.GetNode(n);
			const Node& parent = node.Up();
			ContractionLocal(Rho, Psi[parent], Chi[parent], node, S_opt);
		}
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const TTBasis& basis) {
		Contraction(Rho, Psi, Chi, basis, &S);
	}

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const TTBasis& basis) {
		MatrixTree<T> Rho(basis);
		Contraction(Rho, Psi, Chi, S, basis);
		return Rho;
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TTBasis& basis, bool orthogonal) {
		if (orthogonal) {
			Contraction(Rho, Psi, Psi, basis);
		} else {
			MatrixTree<T> S = DotProduct(Psi, Psi, basis);
			Contraction(Rho, Psi, Psi, S, basis);
		}
	}

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const TTBasis& basis, bool orthogonal) {
		MatrixTree<T> Rho(basis);
		Contraction(Rho, Psi, basis, orthogonal);
		return Rho;
	}
}

#endif //MATRIXTREE_IMPLEMENTATION_H
