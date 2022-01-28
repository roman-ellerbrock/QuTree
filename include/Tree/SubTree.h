//
// Created by Roman Ellerbrock on 1/11/22.
//

#ifndef SUBTREE_H
#define SUBTREE_H
#include "Tree.h"

bool contains(const vector<const Node*>& vec, const Node* probe) {
	for (const Node* node : vec) {
		if (*node == *probe) { return true; }
	}
	return false;
}

vector<const Node*> gatherNodes(const Tree& tree, const vector<size_t>& idx) {
	vector<const Node*> tmp;
	for (size_t i : idx) {
		/// for every node, collect all ancestors
		const Leaf& leaf = tree.leafArray()[i];
		const Node* node = leaf.parent_;
		while(true) {
			if (contains(tmp, node)) { break; }
			tmp.push_back(node);
			if (node->isToplayer()) { break; }
			node = node->parent_;
		}
	}
	return tmp;
}

/**
 * \class SubTree is a mask for nodes in a tree that allows to select a subset of nodes
 */
class SubTree {

public:
	SubTree() = default;

	SubTree(const Tree& tree) {
		initialize(tree);
	}
	~SubTree() = default;

	void initialize(const Tree& tree) {
		nodes_.clear();
		for (Node& node : tree) {
			nodes_.push_back(&node);
		}

		fillEdges(tree);
	}

	void removeTail() {

	}

	void invert() {

	}

	SubTree(const Tree& tree, const vector<size_t>& idx) {
		if (idx.empty()) {
			initialize(tree);
			return;
		}

		/// gather nodes in wrong order
		auto tmp = gatherNodes(tree, idx);

		/// create bottom-up
		for (const Node& node : tree) {
			if (contains(tmp, &node)) {
				nodes_.push_back(&node);
			}
		}

		/// ! idea to avoid n^2 scaling: !
		/// don't go leaf->root but root->leaf (i.e. invert)
		/// then whole vector is top-down.
		/// reverse vector at the very end.

		/// fill edges
		fillEdges(tree);
	}

	void print() const {
		for (const Node* node : nodes_) {
			node->info();
		}
	}

	void fillEdges(const Tree& tree) {
		edges_.clear();
		for (const Edge& edge : tree.edges()) {
			if (contains(nodes_, &edge.from()) && contains(nodes_, &edge.to())) {
				edges_.push_back(&edge);
			}
		}
	}

	vector<const Node*> nodes_;
	vector<const Edge*> edges_;
};


#endif //SUBTREE_H
