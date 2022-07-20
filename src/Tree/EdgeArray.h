//
// Created by Roman Ellerbrock on 11/28/21.
//

#ifndef EDGEARRAY_H
#define EDGEARRAY_H
#include "Edge.h"
#include "NodeArray.h"

class EdgeArray : public vector<Edge> {
public:
	EdgeArray() = default;
	~EdgeArray() = default;

	explicit EdgeArray(const NodeArray& nodeArray) {
		down_.clear();
		for (const Node& node : nodeArray) {
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				down_.emplace_back(Edge(parent, node));
			}
		}

		up_.clear();
		for (const Node& node : nodeArray.reverse_) {
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				up_.emplace_back(Edge(node, parent));
			}
		}

		clear();
		for (const Edge& edge : up_) {
			emplace_back(edge);
		}
		for (const Edge& edge : down_) {
			emplace_back(edge);
		}
	}

	vector<Edge> up_;
	vector<Edge> down_;
};


#endif //EDGEARRAY_H