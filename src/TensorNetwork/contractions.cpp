//
// Created by Roman Ellerbrock on 12/4/21.
//

#include "TensorNetwork/contractions.h"
#include "TensorNetwork/TensorTreeFactory.h"

typedef complex<double> cd;
typedef double d;

/**
 * Node targeting routines
 */

template<typename T>
void apply(Tensor<T>& Ket, const TensorTree<T>& pmat,
	const ProductOperator<T>& P, const Node* node) {

	Tensor<T> tmp(Ket.shape_);
	for (const Edge& edge : pmat.incomingEdges(node)) {
		tmp = matrixTensor(pmat[edge], Ket, edge);
		swap(tmp, Ket);
	}

	for (const Leaf& leaf : node->leaves_) {
		P.apply(tmp, Ket, leaf);
		std::swap(tmp, Ket);
	}
}

template void apply(Tensor<cd>& Ket, const TensorTree<cd>& pmat,
	const ProductOperator<cd>& P, const Node* node);
template void apply(Tensor<d>& Ket, const TensorTree<d>& pmat,
	const ProductOperator<d>& P, const Node* node);


template<typename T>
void apply(Tensor<T>& Ket, const vector<TensorTree<T>>& S,
	const SOP<T>& H, const Node *node) {
	Tensor<T> HKet(Ket.shape_);
	for (size_t l = 0; l < H.size(); ++l) {
		Tensor<T> sKet = Ket;
		apply(sKet, S[l], H[l], node);
		HKet += H.coeff(l) * sKet;
	}
	Ket = HKet;
}

template void apply(Tensor<cd>& Ket,
	const vector<TensorTree<cd>>& S, const SOP<cd>& H,
	const Node *node);
template void apply(Tensor<d>& Ket,
	const vector<TensorTree<d>>& S, const SOP<d>& H,
	const Node *node);


/**
 * Edge targeting routines
 */

template<typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& Ket, const Edge& edge) {
	return matrixTensor(h, Ket, edge.inIdx());
}

template Tensor<cd> matrixTensor(const Tensor<cd>& h, const Tensor<cd>& Ket, const Edge& edge);
template Tensor<d> matrixTensor(const Tensor<d>& h, const Tensor<d>& Ket, const Edge& edge);


template<typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket, const Edge& edge) {
	return contraction(bra, ket, edge.outIdx());
}

template Tensor<cd> contraction(const Tensor<cd>& bra, const Tensor<cd>& ket, const Edge& edge);
template Tensor<d> contraction(const Tensor<d>& bra, const Tensor<d>& ket, const Edge& edge);


template<typename T>
void contraction(TensorTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket,
	const ProductOperator<T>& P, const Edge* edge) {
	apply(Ket, S, P, edge);
	S[edge] = contraction(Bra, Ket, *edge);
}

template void contraction(TensorTree<cd>& S, const Tensor<cd>& Bra, Tensor<cd> Ket,
	const ProductOperator<cd>& P, const Edge* edge);
template void contraction(TensorTree<d>& S, const Tensor<d>& Bra, Tensor<d> Ket,
	const ProductOperator<d>& P, const Edge* edge);


template<typename T>
void contraction(vector<TensorTree<T>>& Svec, const Tensor<T>& Bra,
	Tensor<T> Ket, const SumOfProductsOperator<T>& H, const Edge* edge) {
	for (size_t l = 0; l < Svec.size(); ++l) {
		contraction(Svec[l], Bra, Ket, H[l], edge);
	}
}

template void contraction(vector<TensorTree<cd>>& Svec, const Tensor<cd>& Bra,
	Tensor<cd> Ket, const SumOfProductsOperator<cd>& H, const Edge* edge);
template void contraction(vector<TensorTree<d>>& Svec, const Tensor<d>& Bra,
	Tensor<d> Ket, const SumOfProductsOperator<d>& H, const Edge* edge);


template<typename T>
void apply(Tensor<T>& Ket, const TensorTree<T>& S,
	const ProductOperator<T>& P, const Edge *edge) {

	Tensor<T> sKet(Ket.shape_);
	for (const Edge& preEdge : S.preEdges(edge)) {
		sKet = matrixTensor(S[preEdge], Ket, preEdge);
		std::swap(sKet, Ket);
	}

	const Node& node = edge->from();
	for (const Leaf& leaf : node.leaves_) {
		P.apply(sKet, Ket, leaf);
		std::swap(sKet, Ket);
	}
}

template void apply(Tensor<cd>& Ket, const TensorTree<cd>& S,
	const ProductOperator<cd>& P, const Edge *edge);
template void apply(Tensor<d>& Ket, const TensorTree<d>& S,
	const ProductOperator<d>& P, const Edge *edge);


template<typename T>
void apply(Tensor<T>& Ket, const vector<TensorTree<T>>& S,
	const SOP<T>& H, const Edge *edge) {
	Tensor<T> HKet(Ket.shape_);
	for (size_t l = 0; l < H.size(); ++l) {
		Tensor<T> sKet = Ket;
		apply(sKet, S[l], H[l], edge);
		HKet += H.coeff(l) * sKet;
	}
	Ket = HKet;
}

template void apply(Tensor<cd>& Ket,
	const vector<TensorTree<cd>>& S, const SOP<cd>& H,
	const Edge *edge);
template void apply(Tensor<d>& Ket,
	const vector<TensorTree<d>>& S, const SOP<d>& H,
	const Edge *edge);

/**
 * TensorTree targeting routines
 */

template<typename T>
void apply(TensorTree<T>& Ket, const TensorTree<T>& pmat,
	const ProductOperator<T>& P) {
	for (const Node* node : pmat.nodes_) {
		apply(Ket[node], pmat, P, node);
	}
}

template void apply(TensorTree<cd>& Ket, const TensorTree<cd>& pmat,
	const ProductOperator<cd>& P);
template void apply(TensorTree<d>& Ket, const TensorTree<d>& pmat,
	const ProductOperator<d>& P);


template<typename T>
void contraction(TensorTree<T>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const ProductOperator<T>& P) {
	for (const Edge *edge : S.edges_) {
		contraction(S, Bra[edge], Ket[edge], P, edge);
	}
}

template void contraction(TensorTree<cd>& S, const TensorTree<cd>& Bra,
	TensorTree<cd> Ket, const ProductOperator<cd>& P);
template void contraction(TensorTree<d>& S, const TensorTree<d>& Bra,
	TensorTree<d> Ket, const ProductOperator<d>& P);


template<typename T>
void contraction(vector<TensorTree<T>>& Svec, const TensorTree<T>& Bra,
	TensorTree<T> Ket, const SumOfProductsOperator<T>& H) {

	assert(Svec.size() == H.size());
	for (size_t l = 0; l < Svec.size(); ++l) {
		contraction(Svec[l], Bra, Ket, H[l]);
	}
}

template void contraction(vector<TensorTree<cd>>& Svec, const TensorTree<cd>& Bra, TensorTree<cd> Ket,
	const SumOfProductsOperator<cd>& H);
template void contraction(vector<TensorTree<d>>& Svec, const TensorTree<d>& Bra, TensorTree<d> Ket,
	const SumOfProductsOperator<d>& H);

template<typename T>
vector<TensorTree<T>> contraction(const TensorTree<T>& Bra,
	TensorTree<T> Ket, const SumOfProductsOperator<T>& H) {
	vector<TensorTree<T>> Svec = matrixTree<T>(Bra, H);
	contraction(Svec, Bra, Ket, H);
	return Svec;
} /// untested, generates dense tree

template vector<TensorTree<cd>> contraction(const TensorTree<cd>& Bra,
	TensorTree<cd> Ket, const SumOfProductsOperator<cd>& H);
template vector<TensorTree<d>> contraction(const TensorTree<d>& Bra,
	TensorTree<d> Ket, const SumOfProductsOperator<d>& H);

/*
template<typename T>
void apply(TensorTree<T>& Psi, const vector<TensorTree<T>>& Hmat,
	const SOP<T>& H) {

	TensorTree<T> HPsi(Psi);
	for (size_t l = 0; l < Hmat.size(); ++l) {
		TensorTree<T> sPsi(Psi);
		apply(sPsi, Hmat[l], H[l]);
		for (const Node* node : HPsi.nodes_) {
			HPsi[node] += H.coeff(l) * sPsi[node];
		}
	}
	Psi = HPsi;
}

template void apply(TensorTree<cd>& Psi, const vector<TensorTree<cd>>& Hmat,
	const SOP<cd>& H);
template void apply(TensorTree<d>& Psi, const vector<TensorTree<d>>& Hmat,
	const SOP<d>& H);
*/

template<typename T>
double residual(const TensorTree<T>& Psi1, const TensorTree<T>& Psi2, const Tree& tree) {
	TensorTree<T> S11 = matrixTree<T>(tree);
	TensorTree<T> S22 = matrixTree<T>(tree);
	TensorTree<T> S12 = matrixTree<T>(tree);
	contraction(S11, Psi1, Psi1);
	contraction(S22, Psi2, Psi2);
	contraction(S12, Psi1, Psi2);
	const Tensor<T>& s11 = S11[tree.root()];
	const Tensor<T>& s22 = S22[tree.root()];
	const Tensor<T>& s12 = S12[tree.root()];
	auto diff = s11 + s22 - s12 - adjoint(s12);
	return sqrt(abs(nrm2(diff)));
}

template double residual(const TensorTree<cd>& Psi1, const TensorTree<cd>& Psi2, const Tree& );
template double residual(const TensorTree<d>& Psi1, const TensorTree<d>& Psi2, const Tree&);

