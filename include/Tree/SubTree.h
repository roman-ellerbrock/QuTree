//
// Created by Roman Ellerbrock on 1/11/22.
//

#ifndef SUBTREE_H
#define SUBTREE_H
#include "Tree.h"

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

	bool contains(const Node* probe) {
		for (const Node* node : *this) {
			if (*node == *probe) { return true; }
		}
		return false;
	}

	SubTree(const Tree& tree, const vector<size_t>& idx) {
		if (idx.empty()) { initialize(tree); }
		for (size_t i : idx) {
			cout << i << endl;
			cout << tree.leafArray().size() << endl;
			const Leaf& leaf = tree.leafArray()[i];
			leaf.info();
			const Node* node = leaf.parent_;
			node->info();
			while(true) {
				node->info();
				if (contains(node)) { break; }
				push_back(node);
				if (node->isToplayer()) { break; }
				node = node->parent_;
			}
		}
	}

	void print() const {
		for (const Node* node : *this) {
			node->info();
		}
	}

};


#endif //SUBTREE_H
