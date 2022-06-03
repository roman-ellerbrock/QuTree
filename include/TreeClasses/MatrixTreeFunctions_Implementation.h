//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef MATRIXTREE_IMPLEMENTATION_H
#define MATRIXTREE_IMPLEMENTATION_H
#include "MatrixTreeFunctions.h"
#include "Core/TensorBLAS.h"

namespace TreeFunctions {
////////////////////////////////////////////////////////////////////////
/// General DotProduct for Tensor Trees
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void dotProductLocal(MatrixTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket, const Node& node) {
		if (!node.isBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				const Node& child = node.child(k);
				Ket = matrixTensorBLAS(S[child], Ket, k);
//				Ket = matrixTensor(S[child], Ket, k);
			}
		}
		size_t last = node.shape().lastIdx();
//		contraction(S[node], Bra, Ket, last);
		contractionBLAS(S[node], Bra, Ket, last);

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

	template<typename T>
	double residual(const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		MatrixTree<T> S12 = dotProduct(Bra, Ket, tree);
		MatrixTree<T> S11 = dotProduct(Bra, Bra, tree);
		MatrixTree<T> S22 = dotProduct(Ket, Ket, tree);
		const Node& top = tree.topNode();
		Matrix<T> r = S11[top] + S22[top] - S12[top] - S12[top].adjoint();
		return r.frobeniusNorm();
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
//					Ket = matrixTensor(S[child], Ket, k);
					Ket = matrixTensorBLAS(S[child], Ket, k);
				}
			}
		}

		//Ket = multStateAB(Rho[parent], Ket);
		Ket = matrixTensorBLAS(Rho[parent], Ket, Ket.shape().lastIdx());

		contractionBLAS(Rho[node], Bra, Ket, child_idx);
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
