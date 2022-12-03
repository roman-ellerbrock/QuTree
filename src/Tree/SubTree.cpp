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

/*vector<const Node*> gatherNodes(const Tree& tree, const vector<size_t>& idx) {
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
}*/

vector<const Node*> gatherNodes(const Tree& tree, const vector<size_t>& idx) {
	vector<const Node*> tmp;
	for (size_t i : idx) {
		/// for every node, collect all ancestors
		const Leaf& leaf = tree.leafArray()[i];
		const Node* node = leaf.parent_;
		/// present leaf -> root
		vector<const Node*> slice;
		while(true) {
			if (contains(tmp, node)) { break; }
			slice.push_back(node);
			if (node->isToplayer()) { break; }
			node = node->parent_;
		}
		/// revert, i.e. make root -> leaf
		std::reverse(slice.begin(), slice.end());
		for (const Node* node : slice) {
			tmp.push_back(node);
		}
	}
	return tmp;
}

void SubTree::fillEdges(const Tree& tree) {
	edges_.clear();
	/// for every node in the mask, add the corresponding edge
	for (const Node* node : nodes_) {
		if (node->isToplayer()) { continue; }
		const Node* parent = node->parent_;
		if (nodes_.contains(parent->address_)) {
			Edge helper(node, parent);
			size_t uid = helper.local_address();
			const Edge& edge = tree.edges().up_[uid];
			//size_t address = helper.address();
			edges_.push_back(edge.address(), &edge);
		}
	}

	vector<Edge> downEdges;
	for (auto it = edges_.rbegin(); it != edges_.rend(); ++it) {
		const Edge* edge = *it;
		Edge downEdge(edge->to(), edge->from());
		downEdges.push_back(downEdge);
	}

	for (Edge& tmp : downEdges) {
		size_t uid = tmp.local_address();
		const Edge& edge = tree.edges().down_[uid];
		edges_.push_back(edge.address(), &edge);
	}
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

ostream& operator<<(ostream& os, const SubTree& stree) {
	os << "Nodes:\n";
	for (const Node* node : stree.nodes_) {
		node->info(os);
	}
	os << "Edges:\n";
	for (const Edge* edge : stree.edges_) {
		edge->info(os);
	}
	return os;
}
