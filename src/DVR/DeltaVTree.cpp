//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "DVR/DeltaVTree.h"

void DeltaVTree::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const TensorShape& shape = node.shape();
			size_t dim = shape.lastDimension();
			TensorShape new_shape({dim, dim, dim, dim});
			attributes_.emplace_back(new_shape);
		}
	}
}

void DeltaVTree::print(const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const auto& deltaV = operator[](node);
			node.info();
			deltaV.print();
		}
	}
}
