//
// Created by Roman Ellerbrock on 4/4/22.
//

#ifndef APPLYSCF_H
#define APPLYSCF_H
#include "contractions.h"

template <typename T>
void applyIteration(TensorTree<T>& HPsi, vector<TensorTree<T>>& Hmat,
	const TensorTree<T>& Psi, const SOP<T>& H);

template <typename T>
void apply(TensorTree<T>& HPsi, vector<TensorTree<T>>& mat,
	const TensorTree<T>& Psi, const SOP<T>& H, size_t n_iter);


#endif //APPLYSCF_H
