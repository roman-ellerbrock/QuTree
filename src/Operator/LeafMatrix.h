//
// Created by Roman Ellerbrock on 2020-01-24.
//
#ifndef LEAFMATRIX_H
#define LEAFMATRIX_H
#include "LeafOperator.h"
#include "Tree/Leaf.h"

Tensorcd toMatrix(const LeafOperatorcd& h, const Leaf& leaf);
Tensord toMatrix(const LeafOperatord& h, const Leaf& leaf);

template<typename T>
class LeafMatrix: public LeafOperator<T>
	/**
	 * \class LeafMatrix
	 * \ingroup Operators
	 * \brief This class allows to create LeafOperators from (factor) Matrices.
	 */
{
public:
	LeafMatrix() = default;

	explicit LeafMatrix(Tensor<T> h, bool adjoint = false);

	~LeafMatrix() = default;

	virtual void apply(const BasisAPI& basis, Tensor<T>& hA,
		const Tensor<T>& A) const override;

	const Tensor<T>& matrix() const { return h_; }

private:
	Tensor<T> h_;
};

typedef LeafMatrix<complex<double>> LeafMatrixcd;

typedef LeafMatrix<double> LeafMatrixd;

#endif //LEAFMATRIX_H
