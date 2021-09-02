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

template <typename T>
void toTensor(Tensor<T>& A, const MLO<T>& M, size_t part, const Leaf& leaf);
template <typename T>
void toTensor(Tensor<T>& A, const Matrix<T>& M, size_t part, const Leaf& leaf);
template <typename T>
Matrix<T> toMatrix(const MLO<T>& M, const Leaf& leaf);
template <typename T>
Matrix<T> toMatrix(const Tensor<T>& B, size_t l, const Leaf& leaf);

template <typename T>
class TensorTreeOperator : public TensorTree<T> {
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
	TensorTreeOperator(const MLO<T>& M, const Tree& optree);
	TensorTreeOperator(const SOP<T>& S, const Tree& optree);

	~TensorTreeOperator() = default;

	void occupy(const Tree& optree);

	void occupy(const Tree& optree, mt19937& gen);

	void print(const Tree& optree) const;

	void setLeafOperator(const Matrix<T>& M,
		size_t operator_idx, const Node& node);

	void setLeafOperator(const LeafMatrix<T>& M,
		size_t operator_idx, const Node& node) {
		this->setLeafOperator(M.matrix(), operator_idx, node);
	}

	using TensorTree<T>::attributes_;
	using TensorTree<T>::operator[];
};

template <typename T>
TensorTreeOperator<T> product(const MLO<T>& M, TensorTreeOperator<T> A, const Tree& tree);

typedef TensorTreeOperator<double> TensorTreeOperatord;
typedef TensorTreeOperator<complex<double>> TensorTreeOperatorcd;

#endif //TENSORTREEOPERATOR_H
