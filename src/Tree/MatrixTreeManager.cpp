//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "MatrixTreeManager.h"


namespace MatrixTreeManager {

	template<typename T>
	Matrix<T> DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const TTBasis& basis) {
		for (const Node& node : basis) {
			CalculateLayer(Psi[node], Chi[node], node);
		}
		const Node& topnode = basis.TopNode();
		return S[topnode];
	}

	template<typename T>
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Phi, Tensor<T> AChi, const Node& node) {
		// Get references to the ACoefficients at each node
		if (!node.IsBottomlayer()) {
			for (int k = 0; k < node.nChildren(); k++) {
				// Get overlap-matrix from down_ under
				const Node& child = node.Down(k);
				const Matrix<T>& spo = S[child];
				// Apply it to the right-hand side
				AChi = spo * AChi;
				AChi = multAB(spo, AChi, k);
			}
		}
		DotProduct(S[node], Phi, AChi);
	}

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const MatrixTree<T>& S, const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const TTBasis& basis) {
		assert(Rho.size() == basis.nNodes());
		assert(S.size() == basis.nNodes());
		assert(Psi.size() == basis.nNodes());
		assert(Chi.size() == basis.nNodes());

		for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
			const Node& node = basis.GetNode(n);
			const Node& parent = node.Up();
			ContractionLocal(Rho, S, Psi[parent], Chi[parent], node);
		}
	}

	template<typename T>
	void ContractionLocal(MatrixTree<T>& Rho, const MatrixTree<T>& S, const Tensor<T>& Bra,
		Tensor<T> Ket, const Node& node) {
		assert(!node.IsToplayer());
		const Node& parent = node.Up();

		// Transform sublying basis
		auto child_idx = (size_t) node.ChildIdx();
		for (size_t k = 0; k < parent.nChildren(); ++k) {
			if (k != child_idx) {
				const Node& child = parent.Down(k);
				Ket = multAB(S[child], Ket, k);
			}
		}

		// Transform with upper hole matrix
		if (!parent.IsToplayer()) {
			Ket = multStateAB(Rho[parent], Ket);
		}

		// Calculate hole-product and save to attributes_
		HoleProduct(Rho[node], Bra, Ket, child_idx);
	}

}

