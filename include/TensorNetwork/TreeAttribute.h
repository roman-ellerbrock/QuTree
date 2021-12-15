//
// Created by Roman Ellerbrock on 12/3/21.
//

#ifndef TREEATTRIBUTE_H
#define TREEATTRIBUTE_H
#include "Tree/Tree.h"


template<class C>
class TreeAttribute: public Tree {
public:
	TreeAttribute() = default;
	~TreeAttribute() = default;

	explicit TreeAttribute(const Tree& tree)
		: Tree(tree),
		  nodes_(tree.nNodes()),
		  upEdges_(tree.nEdges()),
		  downEdges_(tree.nEdges()) {
	}

	const C& operator[](const Node& node) const {
		return nodes_[node.address_];
	}

	C& operator[](const Node& node) {
		return nodes_[node.address_];
	}

	const C& operator[](const Edge& edge) const {
		if (edge.isUpEdge()) {
			return upEdges_[edge.address()];
		} else {
			return downEdges_[edge.address()];
		}
	}

	C& operator[](const Edge& edge) {
		if (edge.isUpEdge()) {
			return upEdges_[edge.address()];
		} else {
			return downEdges_[edge.address()];
		}
	}

	vector<C> nodes_;
	vector<C> upEdges_;
	vector<C> downEdges_;
};


#endif //TREEATTRIBUTE_H
