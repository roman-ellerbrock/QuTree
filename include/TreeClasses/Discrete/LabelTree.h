//
// Created by Roman Ellerbrock on 1/23/23.
//

#ifndef LABELTREE_H
#define LABELTREE_H
#include "TreeShape/Tree.h"
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/Discrete/SymmetricSCF.h"
#include "TreeClasses/Discrete/U1Symmetry.h"

class LabelTree {
public:
	LabelTree() = default;
	~LabelTree() = default;

	LabelTree(const Tree& tree, const vector<Labels>& leaf_labels, const Range& range) {
		initialize(tree, leaf_labels, range);
	}

	void initialize(const Tree& tree, const vector<Labels>& leaf_labels, const Range& range) {
		up_.clear();
		for (const Node& node: tree) {
			up_.emplace_back(Labels());
			down_.emplace_back(Labels());
		}

		for (const Node& node: tree) {
			if (node.isBottomlayer()) {
				up_[node] = leaf_labels[node.getLeaf().mode()];
			} else if (!node.isToplayer()) {
				const Node& hole = node.parent();
				up_[node] = label_combinations(up_, down_, range, node, hole);
			}
		}

		down_.clear();
		for (int i = tree.nNodes() - 1; i >= 0; --i) {
			const Node& hole = tree.getNode(i);
			if (hole.isToplayer()) { continue; }
			const Node& node = hole.parent();
			down_[hole] = label_combinations(up_, down_, range, node, hole);
		}
	}

	void print(const Tree& tree) {
		cout << "up:\n";
		for (const Node& node: tree) {
			node.info();
			cout << up_[node] << endl;
		}

		cout << "down:\n";
		for (const Node& node: tree) {
			if (node.isToplayer()) { continue; }
			node.info();
			cout << down_[node] << endl;
		}
	}

	NodeAttribute<Labels> up_;
	NodeAttribute<Labels> down_;
};

using BlockTensor = map<Label, ConfigurationTensor<>>; /// label -> tensor

using LabelDimension = map<Label, size_t>;

ostream& operator<<(ostream& os, const LabelDimension& dim);

size_t labelDimension(const Label& label, size_t dimension, size_t n_hole);

class LabelDimensionTree {
public:
	LabelDimensionTree() = default;
	~LabelDimensionTree() = default;

	LabelDimensionTree(const LabelTree& label_tree, size_t max_dim, const Tree& tree) {
		initialize(label_tree, max_dim, tree);
	}

	LabelDimension createDimension(const NodeAttribute<Labels>& labels,
		size_t nbox, size_t max_dim, const Node& node) {
		map<Label, size_t> dim;
		for (const Label& label: labels[node]) {
			dim[label] = labelDimension(label, max_dim, nbox);
		}
		return dim;
	}

	void initialize(const LabelTree& label_tree, size_t max_dim, const Tree& tree) {
		up_.attributes_.resize(tree.nNodes());
		down_.attributes_.resize(tree.nNodes());

		for (const Node& node: tree) {
			size_t nbox = node.nLeaves();
			up_[node] = createDimension(label_tree.up_, nbox, max_dim, node);
		}

		size_t ntot = tree.topNode().nLeaves();
		for (int n = tree.nNodes() - 1; n >= 0; --n) {
			const Node& node = tree.getNode(n);
			size_t nbox = ntot - node.nLeaves();
			down_[node] = createDimension(label_tree.down_, nbox, max_dim, node);
		}
	}

	void print(const Tree& tree) const {
		cout << "up:\n";
		for (const Node& node: tree) {
			node.info();
			cout << up_[node] << endl;
		}

		cout << "down:\n";
		for (int n = tree.nNodes() - 1; n >= 0; --n) {
			const Node& node = tree.getNode(n);
			node.info();
			cout << down_[node] << endl;
		}
	}

	NodeAttribute<LabelDimension> up_;
	NodeAttribute<LabelDimension> down_;
};

TensorShape LabelShape(const TensorShape& shape, const Label& label);

#endif //LABELTREE_H
