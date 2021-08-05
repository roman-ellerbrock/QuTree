//
// Created by Roman Ellerbrock on 5/23/20.
//

#ifndef TENSOROPERATORTREE_H
#define TENSOROPERATORTREE_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/MultiLeafOperator.h"
#include "TreeOperators/SumOfProductsOperator.h"

class TensorOperatorTree : public TensorTreed {
public:
	TensorOperatorTree() = default;
	explicit TensorOperatorTree(const Tree& tree);
	TensorOperatorTree(const MLOd& M, const Tree& tree);
	TensorOperatorTree(const SOPd& S, const Tree& tree);

	~TensorOperatorTree() = default;

	void occupy(const Tree& tree);

	void print(const Tree& tree) const;

	void setLeafOperator(const Matrixd& M,
		size_t operator_idx, const Node& node);

	void setLeafOperator(const LeafMatrixd& M,
		size_t operator_idx, const Node& node) {
		setLeafOperator(M.matrix(), operator_idx, node);
	}

};


#endif //TENSOROPERATORTREE_H
