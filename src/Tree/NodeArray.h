//
// Created by Roman Ellerbrock on 11/28/21.
//

#ifndef NODEARRAY_H
#define NODEARRAY_H
#include "Tree/Node.h"

class NodeArray: public vector<reference_wrapper<Node>> {
public:
	NodeArray() = default;

	explicit NodeArray(Node& root) {
		clear();
		Node *node = &root;
		while (node) {
			emplace_back(*node);
			node = sweep(node);
		}
		for (size_t i = 0; i < size(); ++i) {
			operator[](i).get().address_ = i;
		}

		reverse_.clear();
		for (auto it = rbegin(); it != rend();++it) {
			reverse_.emplace_back(*it);
		}
	}

	vector<reference_wrapper<Node>> reverse_;
};



#endif //NODEARRAY_H
