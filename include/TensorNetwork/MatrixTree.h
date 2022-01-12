//
// Created by Roman Ellerbrock on 12/24/21.
//

#ifndef MATRIXTREE_H
#define MATRIXTREE_H
#include "TensorTree.h"
#include "Operator/ProductOperator.h"

template <typename T>
class MatrixTree : public TensorTree<T> {
public:
	MatrixTree() = default;
	~MatrixTree() = default;

	ProductOperator<T> p_;
};

#endif //MATRIXTREE_H
