//
// Created by Roman Ellerbrock on 2/10/20.
//

#ifndef TREEIO_H
#define TREEIO_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/TensorTree.h"

namespace IOTree {

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os = cout);

	template <typename T>
	void Leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os = cout);

}

#endif //TREEIO_H