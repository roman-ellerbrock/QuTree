//
// Created by Roman Ellerbrock on 4/4/22.
//

#include "ApplySCF.h"
#include "contractions.h"

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
