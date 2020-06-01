//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef TREEFACTORY_H
#define TREEFACTORY_H
#include "TreeShape/Tree.h"

namespace TreeFactory {
	Tree BalancedTree(size_t num_leaves, size_t dim_leaves, size_t dim_nodes);
	Tree UnbalancedTree(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType);
	Tree OperatorTree(const Tree& tree);

}

#endif //TREEFACTORY_H
