//
// Created by Roman Ellerbrock on 11/28/21.
//

#ifndef EDGEARRAY_H
#define EDGEARRAY_H
#include "Edge.h"
#include "NodeArray.h"

class EdgeArray: public vector<Edge> {
public:
	EdgeArray() = default;
	~EdgeArray() = default;

	explicit EdgeArray(const NodeArray& nodeArray) {
		clear();
		for (const Node& node : nodeArray) {
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				emplace_back(Edge(node, parent));
			}
		}
	}
};


#endif //EDGEARRAY_H
