//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef TREEFACTORY_H
#define TREEFACTORY_H
#include "TreeShape/Tree.h"

namespace TreeFactory {
	vector<Node> partition(const vector<Node>& nodes, size_t n_partition, size_t dim_node);
	vector<Node> bottomlayerNodes(size_t num_leaves, size_t dim_leaves, size_t dim_nodes, size_t leaf_type);
	Tree balancedTree(size_t num_leaves, size_t dim_leaves = 2,
		size_t dim_nodes = 2, size_t dim_inc = 0, size_t leaf_type = 0,
		double omega = 1., double r0 = 0., double wfr0 = 1., double wfomega = 0.);
	Tree unbalancedTree(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType);
	Tree operatorTree(const Tree& tree, size_t bottom = 0);
	map<size_t, size_t> leaves_staggered_integers(size_t num_integer, size_t num_bits);

}

#endif //TREEFACTORY_H
