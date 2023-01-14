//
// Created by Roman Ellerbrock on 1/5/23.
//

#ifndef BLOCKTREESHAPE_H
#define BLOCKTREESHAPE_H
#include "TreeShape/Tree.h"
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/Discrete/SymmetricSCF.h"

using BlockTensorShape = vector<TensorShape>;
using BlockTreeShape = NodeAttribute<BlockTensorShape>;

/// Symmetry label, e.g. particle number
using Label = size_t;
using Labels = Configuration<Label>;
using LabelTree = NodeAttribute<Labels>;

using BlockTensor = map<size_t, ConfigurationTensor<>>; /// label -> tensor

/// range of allowed labels, e.g. particle number <= 5 and >= 0
class Range {
public:
	Range(Label min, Label max) : min_(min), max_(max) {}
	Label min_{0};
	Label max_{1};

	bool isAllowed(size_t x) const {
		return ((x >= min_) && (x <= max_));
	}
};

Labels combine(const Labels& L, const Labels& R, const Range& range);

void add(Label& label, const Label& summand, const Range& range);
Labels label_combinations(const LabelTree& up, const LabelTree& down,
	const Range& range, const Node& node, const Node& hole);

class BlockTree {
public:
	BlockTree() = default;
	~BlockTree() = default;
	BlockTree(const Tree& tree) { initialize(tree); }

	void initialize(const Tree& tree) {
		labels_up_.clear();
		labels_down_.clear();
		for (const Node& node : tree) {
			labels_up_.emplace_back(Labels());
			labels_down_.emplace_back(Labels());
		}
	}

	void initLabels(const Tree& tree, const vector<Labels>& leaf_labels, const Range& range) {
		for (const Node& node: tree) {
			if (node.isBottomlayer()) {
				labels_up_[node] = leaf_labels[node.getLeaf().mode()];
			} else if (!node.isToplayer()) {
				const Node& hole = node.parent();
				labels_up_[node] = label_combinations(labels_up_, labels_down_, range, node, hole);
			}
		}

		for (int i = tree.nNodes() - 1; i >= 0; --i) {
			const Node& hole = tree.getNode(i);
			if (hole.isToplayer()) { continue; }
			const Node& node = hole.parent();
			labels_down_[hole] = label_combinations(labels_up_, labels_down_, range, node, hole);
		}
	}

	void InitShapes(const Tree& tree);

	void print(const Tree& tree) {
		cout << "up:\n";
		for (const Node& node : tree) {
			node.info();
			cout << labels_up_[node] << endl;
		}

		cout << "down:\n";
		for (const Node& node : tree) {
			if (node.isToplayer()) { continue; }
			node.info();
			cout << labels_down_[node] << endl;
		}
	}

	LabelTree labels_up_;
	LabelTree labels_down_;
	BlockTreeShape shapes_;
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