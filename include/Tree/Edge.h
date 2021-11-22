//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef EDGE_H
#define EDGE_H
#include "TreeShape/Node.h"


class Edge {
	const Node* down_;
	const Node* up_;

public:
	Edge(const Node& down, const Node& up) : down_(&down), up_(&up) {}

	const Node& down() const { return *down_; }
	const Node& up() const { return *up_; }

	size_t downIdx() const { return down().nChildren(); }
	size_t upIdx() const { return down().childIdx(); }

};


#endif //EDGE_H
