//
// Created by Roman Ellerbrock on 1/11/22.
//

#include "Tree/SubTree.h"

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

ostream& operator<<(ostream& os, const SubTreeParameters& par) {
	os << "Leaves: ";
	for (auto x : par.leaves_) {
		os << x << ", ";
	}
	os << endl;

	cout << "Nodes active: " << par.activeNodes_ << endl;
	cout << "Edges active: " << par.activeEdges_ << endl;

	return os;
}
