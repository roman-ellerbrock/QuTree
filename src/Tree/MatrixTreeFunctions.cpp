//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "MatrixTreeFunctions.h"


namespace MatrixTreeFunctions {

	template<typename T>
	Matrix<T> DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const TTBasis& basis) {
		for (const Node& node : basis) {
			DotProductLocal(S, Psi[node], Chi[node], node);
		}
		const Node& topnode = basis.TopNode();
		return S[topnode];
	}

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
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const TTBasis& basis) {
		assert(Rho.size() == basis.nNodes());
		assert(S.size() == basis.nNodes());
		assert(Psi.size() == basis.nNodes());
		assert(Chi.size() == basis.nNodes());

		for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
			const Node& node = basis.GetNode(n);
			const Node& parent = node.Up();
			ContractionLocal(Rho, Psi[parent], Chi[parent], S, node);
		}
	}

	template<typename T>
	void ContractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
		const MatrixTree<T>& S, const Node& node) {
		assert(!node.IsToplayer());

		const Node& parent = node.Up();
		auto child_idx = (size_t) node.ChildIdx();
		for (size_t k = 0; k < parent.nChildren(); ++k) {
			if (k != child_idx) {
				const Node& child = parent.Down(k);
				Ket = multAB(S[child], Ket, k);
			}
		}

		if (!parent.IsToplayer()) {
			Ket = multStateAB(Rho[parent], Ket);
		}

		mHoleProduct(Rho[node], Bra, Ket, child_idx);
	}

	typedef complex<double> cd;

	template void DotProductLocal(MatrixTree<cd>& S, const Tensor<cd>& Phi, Tensor<cd> AChi, const Node& node);
	template Matrix<cd>
	DotProduct(MatrixTree<cd>& S, const TensorTree<cd>& Psi, const TensorTree<cd>& Chi, const TTBasis& basis);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const TensorTree<cd>& Chi,
		const MatrixTree<cd>& S, const TTBasis& basis);
	template void
	ContractionLocal(MatrixTree<cd>& Rho, const Tensor<cd>& Bra, Tensor<cd> Ket, const MatrixTree<cd>& S, const Node& node);

	typedef double d;

	template void DotProductLocal<d>(MatrixTree<d>& S, const Tensor<d>& Phi, Tensor<d> AChi, const Node& node);
	template Matrix<d>
	DotProduct<d>(MatrixTree<d>& S, const TensorTree<d>& Psi, const TensorTree<d>& Chi, const TTBasis& basis);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const TensorTree<d>& Chi,
		const MatrixTree<d>& S, const TTBasis& basis);
	template void ContractionLocal(MatrixTree<d>& Rho, const Tensor<d>& Bra, Tensor<d> Ket,
		const MatrixTree<d>& S, const Node& node);
}
