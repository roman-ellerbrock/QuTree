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
	explicit TensorOperatorTree(const Tree& optree);
	TensorOperatorTree(const Tree& optree, mt19937& gen)
	: TensorOperatorTree(optree) {
		occupy(optree, gen);
	}
	TensorOperatorTree(const MLOd& M, const Tree& optree);
	TensorOperatorTree(const SOPd& S, const Tree& optree);

	~TensorOperatorTree() = default;

	void occupy(const Tree& optree);

	void occupy(const Tree& optree, mt19937& gen);

	void print(const Tree& optree) const;

	void setLeafOperator(const Matrixd& M,
		size_t operator_idx, const Node& node);

	void setLeafOperator(const LeafMatrixd& M,
		size_t operator_idx, const Node& node) {
		setLeafOperator(M.matrix(), operator_idx, node);
	}

};


#endif //TENSOROPERATORTREE_H
