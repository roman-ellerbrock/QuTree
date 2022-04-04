//
// Created by Roman Ellerbrock on 1/11/22.
//

#ifndef SUBTREE_H
#define SUBTREE_H
#include "Tree.h"
#include "sparse_vector.h"

bool contains(const vector<const Node *>& vec, const Node *probe);

vector<const Node *> gatherNodes(const Tree& tree, const vector<size_t>& idx);


class SubTree;
vector<const Node *> gatherNodes(const SubTree& tree, const vector<size_t>& idx);

struct SubTreeParameters {
	SubTreeParameters(const vector<size_t>& leaves = {},
		bool activeNodes = true, bool activeEdges = true)
		: leaves_(leaves), activeNodes_(activeNodes), activeEdges_(activeEdges) {}

	~SubTreeParameters() = default;

	vector<size_t> leaves_ = {};
	bool activeNodes_ = true;
	bool activeEdges_ = true;
};

ostream& operator<<(ostream& os, const SubTreeParameters& par);

/**
 * \class SubTree is a mask for nodes in a tree that allows to select a subset of nodes
 */
class SubTree {

public:
	SubTree() = default;
	~SubTree() = default;

	void removeTail() {
		// @TODO: add in the future
		exit(12);
	}

	void invert() {
		// @TODO: add in the future
		exit(12);
	}

	explicit SubTree(const Tree& tree, SubTreeParameters par = SubTreeParameters()) {
		if (par.leaves_.empty()) {
			for (size_t i = 0; i < tree.leafArray().size(); ++i) {
				par.leaves_.push_back(i);
			}
		}

		/// gather nodes in wrong order
		auto tmp = gatherNodes(tree, par.leaves_);

		/// create bottom-up
		for (const Node& node : tree) {
			if (contains(tmp, &node)) {
				nodes_.push_back(node.address_, &node);
			}
		}

		///
		for (auto i : par.leaves_) {
			leaves_[i] = &(tree.leafArray()[i]);
		}

		/// ! idea to avoid n^2 scaling: !
		/// don't go leaf->root but root->leaf (i.e. invert)
		/// then whole vector is top-down.
		/// reverse vector at the very end.

		/// fill edges
		if (par.activeEdges_) { fillEdges(tree); }

		if (!par.activeNodes_) { nodes_.clear(); }
	}

/*	SubTree(const SubTree& tree, SubTreeParameters par = SubTreeParameters()) {
		if (par.leaves_.empty()) {
			for (size_t i = 0; i < tree.leaves_.size(); ++i) {
				par.leaves_.push_back(i);
			}
		}

		/// gather nodes in wrong order
		auto tmp = gatherNodes(tree, par.leaves_);

		/// create bottom-up
		for (const Node* node : tree.nodes_) {
			if (contains(tmp, node)) {
				nodes_.push_back(node->address_, node);
			}
		}

		///
//		for (pair<size_t, const Leaf*> iPair : par.leaves_) {
		for (size_t i : par.leaves_) {
//			size_t i = iPair.first;
			leaves_[i] = tree.leaves_.at(i);
		}

		/// ! idea to avoid n^2 scaling: !
		/// don't go leaf->root but root->leaf (i.e. invert)
		/// then whole vector is top-down.
		/// reverse vector at the very end.

		/// fill edges
		if (par.activeEdges_) { fillEdges(tree); }

		if (!par.activeNodes_) { nodes_.clear(); }

	}*/

	virtual void print() const {
		for (const Node *node : nodes_) {
			node->info();
		}
	}

	void fillEdges(const Tree& tree) {
		edges_.clear();
		for (const Edge& edge : tree.edges()) {
			if (nodes_.contains(edge.from().address_) && nodes_.contains(edge.to().address_)) {
				edges_.push_back(edge.address(), &edge);
			}
		}
	}

	[[nodiscard]] vector<Edge> preEdges(const Edge *edge) const {
		/// Get all pre-Edges and filter for the ones in the SubTree
		vector<Edge> preEd = ::preEdges(*edge);
		vector<Edge> pres;
		for (const Edge& e : preEd) {
			if (edges_.contains(e.address())) {
				pres.push_back(e);
			}
		}
		return pres;
	}

	sparse_vector<const Node *> nodes_;
	sparse_vector<const Edge *> edges_;
	map<size_t, const Leaf *> leaves_;
};


#endif //SUBTREE_H
