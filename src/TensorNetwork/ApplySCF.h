//
// Created by Roman Ellerbrock on 4/4/22.
//

#ifndef APPLYSCF_H
#define APPLYSCF_H
#include "contractions.h"

template <typename T>
void applyIteration(TensorTree<T>& HPsi, vector<TensorTree<T>>& Hmat,
	const TensorTree<T>& Psi, const SOP<T>& H);



#endif //APPLYSCF_H
