//
// Created by Roman Ellerbrock on 4/4/22.
//

#include "ApplySCF.h"
#include "contractions.h"

typedef complex<double> cd;
typedef double d;

template <typename T>
void applyIteration(TensorTree<T>& HPsi, vector<TensorTree<T>>& Hmat,
	TensorTree<T>& sPsi, const TensorTree<T>& Psi, const SOP<T>& H) {

	for (const Edge* edge : HPsi.edges_) {
		/// apply to edge
		apply(HPsi, sPsi, Hmat, H, edge);

		// normalize

		/// represent operator at node
//		contraction()
	}

}

template void applyIteration(TensorTree<cd>& HPsi, vector<TensorTree<cd>>& Hmat,
	TensorTree<cd>& sPsi, const TensorTree<cd>& Psi, const SOP<cd>& H);
template void applyIteration(TensorTree<d>& HPsi, vector<TensorTree<d>>& Hmat,
	TensorTree<d>& sPsi, const TensorTree<d>& Psi, const SOP<d>& H);

template <typename T>
void apply(TensorTree<T>& HPsi, vector<TensorTree<T>>& Hmat,
	const TensorTree<T>& Psi, const SOP<T>& H) {

	size_t n_iter = 10;
	for (size_t i = 0; i < n_iter; ++i) {
		contraction(Hmat, Psi, Psi, H);
		apply(Psi, Hmat, H);
		Psi.normalize();
	}

}
