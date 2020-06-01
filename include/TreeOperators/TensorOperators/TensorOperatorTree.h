//
// Created by Roman Ellerbrock on 5/23/20.
//

#ifndef TENSOROPERATORTREE_H
#define TENSOROPERATORTREE_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/LeafMatrix.h"

class TensorOperatorTree : public TensorTreecd {
public:
	TensorOperatorTree() = default;
	explicit TensorOperatorTree(const Tree& tree);
	~TensorOperatorTree() = default;

	void Occupy(const Tree& tree);

	void print(const Tree& tree) const;

	void setLeafOperator(const LeafMatrixcd& M, size_t operator_idx, const Node& node);

private:
};


#endif //TENSOROPERATORTREE_H
