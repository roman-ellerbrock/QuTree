//
// Created by Roman Ellerbrock on 1/5/23.
//

#ifndef BLOCKTREESHAPE_H
#define BLOCKTREESHAPE_H
#include "TreeClasses/Discrete/LabelTree.h"

using BlockTensorShape = vector<TensorShape>;
using BlockTreeShape = NodeAttribute<BlockTensorShape>;

class BlockTree {
public:
	BlockTree() = default;
	~BlockTree() = default;

	explicit BlockTree(const Tree& tree, const Range& range,
		const vector<Labels>& leaf_labels, size_t max_dim) {
		initialize(tree, range, leaf_labels, max_dim);
	}

	void initialize(const Tree& tree, const Range& range,
		const vector<Labels>& leaf_labels, size_t max_dim) {
		labels_.initialize(tree, leaf_labels, range);
		dim_.initialize(labels_, max_dim, tree);
	}

	BlockTensor bottomBlockTensor(const Labels& labels) {
		BlockTensor A;
		size_t i = 0;
		for (const Label& label: labels) {
			ConfigurationTensor<> c = {{i++}};
			A[label] = c;
		}
		return A;
	}

	BlockTensor upperBlockTensor(const Labels& labels, const Node& node) {
		BlockTensor A;
		for (const Label& label: labels) {

		}
	}

	void print(const Tree& tree) {
		labels_.print(tree);
	}

	LabelTree labels_;
	LabelDimensionTree dim_;
	NodeAttribute<BlockTensor> up_;
	NodeAttribute<BlockTensor> down_;
};

#endif //BLOCKTREESHAPE_H

/**
 * stuff todo:
 * - == for node
 * - rename SparseTree to TreeMask
 * - use SparseTree for dense as well
 * - every NodeAttribute has a SparseTree built in
 * 		- enables operator<<
 * -
 *
*/