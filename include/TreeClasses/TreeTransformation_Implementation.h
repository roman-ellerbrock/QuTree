//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef TREETRANSFORMATION_IMPLEMENTATION_H
#define TREETRANSFORMATION_IMPLEMENTATION_H

#include "TreeClasses/TreeTransformation.h"
#include "TreeClasses/SpectralDecompositionTree.h"


namespace TreeFunctions {


	template <typename T>
	TensorTree<T> DotProductNormalization(TensorTree<T> Psi, const Tree& tree) {

		return Psi;
	}

	template <typename T>
	TensorTree<T> DirectionalInvarientRepresentation(TensorTree<T> Psi, const Tree& tree, bool orthogonal) {
		/// Calculate contraction matrixtree
		assert(orthogonal);
		MatrixTree<T> rho = contraction(Psi, tree, orthogonal);
		auto B = sqrt(rho, tree);


		for (const Edge& e : tree.edges()) {
			TransformEdgeDown(Psi, Psi, B[e], e);
		}
	}

	template <typename T>
	void TransformEdgeDown(TensorTree<T>& Chi, const TensorTree<T>& Psi, const Matrix<T>& M, const Edge& e) {
		const Node& node = e.down();
		Chi[node] = tensorMatrix(Psi[node], M, e.downIdx());
	}

	template <typename T>
	TensorTree<T> contractionNormalization(TensorTree<T> Psi, const Tree& tree, bool orthogonal) {
		assert(orthogonal);
		/// transform to Directional-Invariant rep
		MatrixTree<T> rho = contraction(Psi, tree, orthogonal);
		auto B = sqrt(rho, tree);
		auto B_inv = inverse(B, tree);

		for (const Edge& e : tree.edges()) {
			TransformEdgeDown(Psi, Psi, B[e], e);
		}

		auto Chi(Psi);

		for (const Edge& e : tree.edges()) {
			const Node& parent = e.up();
			Chi[e] = matrixTensor(B_inv[e], Psi[parent], e.upIdx());
		}

		return Chi;
	}

	template <typename T>
	void TransformEdgeUp(TensorTree<T>& Chi, const TensorTree<T>& Psi, const Matrix<T>& Mi, const Edge& e) {
		const Node& parent = e.up();
		Chi[parent] = matrixTensor(Mi, Psi[parent], e.upIdx());
	}

	template <typename T>
	void TransformEdge(TensorTree<T>& Chi, const TensorTree<T>& Psi, const Matrix<T>& M,
		const Matrix<T>& M_inv, const Edge& e) {

		const Node& node = e.down();
		Chi[node] = tensorMatrix(Psi[node], M, node.parentIdx());

		const Node& parent = e.up();
		Chi[parent] = matrixTensor(M_inv, Psi[parent], node.childIdx());
	}

	template <typename T>
	void transform(TensorTree<T>& Chi, const TensorTree<T>& Psi, const MatrixTree<T>& M, const MatrixTree<T>& M_inv,
		const Tree& tree) {

		for (const Edge& e : tree.edges()) {
			TransformEdgeUp(Chi, Psi, M_inv[e], e);
			const Node& node = e.down();
			auto x = contraction(Chi[node], Chi[node], node.childIdx());
			node.info();
			x.print();
		}
		getchar();

	}

}

#endif //TREETRANSFORMATION_IMPLEMENTATION_H
