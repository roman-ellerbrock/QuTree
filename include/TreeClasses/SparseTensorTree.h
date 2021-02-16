//
// Created by Roman Ellerbrock on 2/16/21.
//

#ifndef SPARSETENSORTREE_H
#define SPARSETENSORTREE_H
#include "TreeClasses/SparseNodeAttribute.h"
#include "TreeClasses/TensorTree.h"

template<typename T>
class SparseTensorTree: public SparseNodeAttribute<Tensor<T>> {
	using SparseNodeAttribute<Tensor<T>>::attributes_;
public:
	using SparseNodeAttribute<Tensor<T>>::sparseTree;
	using SparseNodeAttribute<Tensor<T>>::operator[];
	using SparseNodeAttribute<Tensor<T>>::initialize;
	using SparseNodeAttribute<Tensor<T>>::size;

	SparseTensorTree(const SparseTree& stree, const Tree& tree)
		: SparseNodeAttribute<Tensor<T>>(stree, tree) {
		initialize(tree);
	}

	SparseTensorTree(shared_ptr<SparseTree>& stree, const Tree& tree)
		: SparseNodeAttribute<Tensor<T>>(stree, tree) {
		initialize(tree);
	}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	SparseTensorTree(const vector<size_t>& leaf_indices, const Tree& tree, bool tail = true, bool inverse_tree = false)
		: SparseNodeAttribute<Tensor<T>>(leaf_indices, tree, tail, inverse_tree) {
		initialize(tree);
	}


	void initialize(const Tree& tree) override;
};

typedef SparseTensorTree<complex<double>> SparseTensorTreecd;

typedef SparseTensorTree<double> SparseTensorTreed;

#endif //SPARSETENSORTREE_H
