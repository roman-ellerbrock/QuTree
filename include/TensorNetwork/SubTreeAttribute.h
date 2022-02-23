//
// Created by Roman Ellerbrock on 2/16/22.
//

#ifndef SUBTREEATTRIBUTE_H
#define SUBTREEATTRIBUTE_H
#include "Tree/SubTree.h"
#include <map>

template<class C>
class SubTreeAttribute: public SubTree {
public:
	SubTreeAttribute() = default;
	~SubTreeAttribute() = default;

	explicit SubTreeAttribute(const Tree& tree, const vector<size_t>& idxs = {})
		: SubTree(tree, idxs) {
	}

	explicit SubTreeAttribute(const SubTree& subtree) : SubTree(subtree) {
	}

	/// Nodes

	const C& operator[](const Node& node) const {
		return nodeAttribute_.at(node.address_);
	}

	const C& operator[](const Node* node) const {
		return nodeAttribute_.at(node->address_);
	}

	C& operator[](const Node& node) {
		return nodeAttribute_[node.address_];
	}

	C& operator[](const Node* node) {
		return nodeAttribute_[node->address_];
	}

	/// Edges

	const C& operator[](const Edge& edge) const {
		return edgeAttribute_.at(edge.address());
	}

	const C& operator[](const Edge* edge) const {
		return edgeAttribute_.at(edge->address());
	}

	C& operator[](const Edge& edge) {
		return edgeAttribute_[edge.address()];
	}

	C& operator[](const Edge* edge) {
		return edgeAttribute_[edge->address()];
	}

	map<size_t, C> nodeAttribute_;
	map<size_t, C> edgeAttribute_;
};



#endif //SUBTREEATTRIBUTE_H
