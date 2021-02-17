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

	~SparseTensorTree() {
		delete[] mem_;
	}

	void initialize(const Tree& tree) override {
		/**
		 * Rationale:
		 * - pool memory allocator
		 */
		attributes_.clear();
		size_t dim = 0;
		for (const Node *const node_ptr : sparseTree()) {
			dim += node_ptr->shape().totalDimension();
		}
		mem_ = (T *) malloc(dim * sizeof(T)); /// allocate memory for all tensors at once
		memset(mem_, 0, dim * sizeof(T));
		T *loc = mem_; /// point to location of next tensor
		for (const Node *const node_ptr : sparseTree()) {
			attributes_.emplace_back(Tensor<T>(node_ptr->shape(), loc, false, false));
			loc += node_ptr->shape().totalDimension();
		}
	}

private:
	T *mem_;
};

typedef SparseTensorTree<complex<double>> SparseTensorTreecd;

typedef SparseTensorTree<double> SparseTensorTreed;

#endif //SPARSETENSORTREE_H
