//
// Created by Roman Ellerbrock on 4/4/22.
//

#include "ApplySCF.h"

typedef complex<double> cd;
typedef double d;

template <typename T>
void applyIteration(TensorTree<T>& HPsi, vector<TensorTree<T>>& mat,
	const TensorTree<T>& Psi, const SOP<T>& H) {

	for (const Edge* edge : HPsi.edges_) {
		/// apply to edge
		const Node& node = edge->from();
		HPsi[edge] = Psi[node];
		apply(HPsi[edge], mat, H, edge);

		/// normalize
		HPsi[edge] = normalize(HPsi[edge], edge);

		/// represent operator at node
		contraction(mat, HPsi[edge], Psi[edge], H, edge);
	}
}

template void applyIteration(TensorTree<cd>& HPsi, vector<TensorTree<cd>>& Hmat,
	const TensorTree<cd>& Psi, const SOP<cd>& H);
template void applyIteration(TensorTree<d>& HPsi, vector<TensorTree<d>>& Hmat,
	const TensorTree<d>& Psi, const SOP<d>& H);


template <typename T>
void apply(TensorTree<T>& HPsi, vector<TensorTree<T>>& mat,
	const TensorTree<T>& Psi, const SOP<T>& H) {

	size_t n_iter = 10;
	for (size_t i = 0; i < n_iter; ++i) {
		applyIteration(HPsi, mat, Psi, H);
	}

}

template void apply(TensorTree<cd>& HPsi, vector<TensorTree<cd>>& mat,
	const TensorTree<cd>& Psi, const SOP<cd>& H);
template void apply(TensorTree<d>& HPsi, vector<TensorTree<d>>& mat,
	const TensorTree<d>& Psi, const SOP<d>& H);
