//
// Created by Roman Ellerbrock on 5/23/20.
//

#ifndef TENSORTREEOPERATOR_H
#define TENSORTREEOPERATOR_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/TensorTree.h"
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/MultiLeafOperator.h"
#include "TreeOperators/SumOfProductsOperator.h"

void toTensor(Tensord& A, const MLOd& M, size_t part, const Leaf& leaf);
void toTensor(Tensord& A, const Matrixd& M, size_t part, const Leaf& leaf);
Matrixd toMatrix(const MLOd& M, const Leaf& leaf);
Matrixd toMatrix(const Tensord& B, size_t l, const Leaf& leaf);

class TensorTreeOperator : public TensorTreed {
	/**
	 * \brief This class is a tree tensor network operator (TTNO)
	 *
	 * \defgroup TTNO
	 *
	 */
public:
	TensorTreeOperator() = default;
	explicit TensorTreeOperator(const Tree& optree);
	TensorTreeOperator(const Tree& optree, mt19937& gen)
	: TensorTreeOperator(optree) {
		occupy(optree, gen);
	}
	TensorTreeOperator(const MLOd& M, const Tree& optree);
	TensorTreeOperator(const SOPd& S, const Tree& optree);

	~TensorTreeOperator() = default;

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

TensorTreeOperator product(const MLOd& M, TensorTreeOperator A, const Tree& tree);

#endif //TENSORTREEOPERATOR_H
