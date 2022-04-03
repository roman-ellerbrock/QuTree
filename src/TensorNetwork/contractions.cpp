//
// Created by Roman Ellerbrock on 12/4/21.
//

#include "TensorNetwork/contractions.h"

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
void apply(TensorTree<T>& Ket, const TensorTree<T>& S, const ProductOperator<T>& P,
	const Edge* edge) {

	Tensor<T> sKet(Ket[edge].shape_);
	for (const Edge& preEdge : S.preEdges(edge)) {
		sKet = matrixTensor(S[preEdge], Ket[edge], preEdge);
		std::swap(sKet, Ket[edge]);
	}

	const Node& node = edge->from();
	for (const Leaf& leaf : node.leaves_) {
		P.apply(sKet, Ket[edge], leaf);
		std::swap(sKet, Ket[edge]);
	}
}

template void apply(TensorTree<cd>& Ket, const TensorTree<cd>& S,
	const ProductOperator<cd>& P, const Edge* edge);
template void apply(TensorTree<d>& Ket, const TensorTree<d>& S,
	const ProductOperator<d>& P, const Edge* edge);

template <typename T>
void contraction(TensorTree<T>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const ProductOperator<T>& P) {
	for (const Edge* edge : S.edges_) {
		apply(Ket, S, P, edge);
		S[edge] = contraction(Bra[edge], Ket[edge], *edge);
	}
}

template void contraction(TensorTree<cd>& S, const TensorTree<cd>& Bra,
	TensorTree<cd> Ket, const ProductOperator<cd>& P);
template void contraction(TensorTree<d>& S, const TensorTree<d>& Bra,
	TensorTree<d> Ket, const ProductOperator<d>& P);

template <typename T>
void contraction(vector<TensorTree<T>>& Svec, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const SumOfProductsOperator<T>& H) {
	assert(Svec.size() == sop.size());
	for (size_t l = 0; l < Svec.size(); ++l) {
		contraction(Svec[l], Bra, Ket, H[l]);
	}
}

template void contraction(vector<TensorTree<cd>>& Svec, const TensorTree<cd>& Bra, TensorTree<cd> Ket,
	const SumOfProductsOperator<cd>& H);
template void contraction(vector<TensorTree<d>>& Svec, const TensorTree<d>& Bra, TensorTree<d> Ket,
	const SumOfProductsOperator<d>& H);

template <typename T>
TensorTree<T> dotProduct(const TensorTree<T>& Bra, TensorTree<T> Ket) {
	TensorTree<T> S;
	contraction(S, Bra, Ket, Bra.shape_.lastIdx());
	return S;
}

template <typename T>
Tensor<T> apply(Tensor<T> Ket, const TensorTree<T>& S, const Node& node) {
	Tensor<T> tmp(Ket.shape_);
	for (const Edge& edge : incomingEdges(node)) {
		tmp = matrixTensor(S[edge], Ket, edge);
		swap(tmp, Ket);
	}
	return Ket;
}

template Tensor<cd> apply(Tensor<cd> Ket, const TensorTree<cd>& S, const Node& );
template Tensor<d> apply(Tensor<d> Ket, const TensorTree<d>& S, const Node&);
