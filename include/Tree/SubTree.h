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

class SubTree : public vector<const Node*> {

public:
	SubTree() = default;

	SubTree(const Tree& tree) {
		initialize(tree);
	}
	~SubTree() = default;

	void initialize(const Tree& tree) {
		clear();
		for (Node& node : tree) {
			push_back(&node);
		}
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
				push_back(&node);
			}
		}

		/// ! idea to avoid n^2 scaling: !
		/// don't go leaf->root but root->leaf (i.e. invert)
		/// then whole vector is top-down.
		/// reverse vector at the very end.
	}

	void print() const {
		for (const Node* node : *this) {
			node->info();
		}
	}

};


#endif //SUBTREE_H
