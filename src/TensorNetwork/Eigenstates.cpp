//
// Created by Roman Ellerbrock on 4/3/22.
//

#include "Eigenstates.h"

template <typename T>
void relaxation(TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree) {

	vector<TensorTree<T>> Hmat = matrixTree(tree, H);

}