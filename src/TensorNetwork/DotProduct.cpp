//
// Created by Roman Ellerbrock on 12/4/21.
//

#include "TensorNetwork/DotProduct.h"

typedef complex<double> cd;
typedef double d;

template <typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& Ket, const Edge& edge) {
	return matrixTensor(h, Ket, edge.inIdx());
}

template Tensor<cd> matrixTensor(const Tensor<cd>& h, const Tensor<cd>& Ket, const Edge& edge);
template Tensor<d> matrixTensor(const Tensor<d>& h, const Tensor<d>& Ket, const Edge& edge);


template <typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket, const Edge& edge) {
	return contraction(bra, ket, edge.outIdx());
}

template Tensor<cd> contraction(const Tensor<cd>& bra, const Tensor<cd>& ket, const Edge& edge);
template Tensor<d> contraction(const Tensor<d>& bra, const Tensor<d>& ket, const Edge& edge);


template <typename T>
void apply(TensorTree<T>& Ket, const TensorTree<T>& S, const Edge& edge) {
	vector<Edge> preedges = preEdges(edge);
	for (const Edge& preEdge : preedges) {
		Ket[edge] = matrixTensor(S[preEdge], Ket[edge], preEdge);
	}
}

template void apply(TensorTree<cd>& Ket, const TensorTree<cd>& S, const Edge& edge);
template void apply(TensorTree<d>& Ket, const TensorTree<d>& S, const Edge& edge);


template <typename T>
void dotProduct(TensorTree<T>& S, const TensorTree<T>& Bra, TensorTree<T> Ket) {
	for (const Edge& edge : S.edges()) {
		apply(Ket, S, edge);
		S[edge] = contraction(Bra[edge], Ket[edge], edge);
	}
}

template void dotProduct(TensorTree<cd>& S, const TensorTree<cd>& Bra, TensorTree<cd> Ket);
template void dotProduct(TensorTree<d>& S, const TensorTree<d>& Bra, TensorTree<d> Ket);
