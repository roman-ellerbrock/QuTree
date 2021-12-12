//
// Created by Roman Ellerbrock on 12/3/21.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "TreeAttribute.h"

template<typename T>
class TensorTree: public TreeAttribute<Tensor<T>> {
public:
	TensorTree() = default;
	~TensorTree() = default;

	explicit TensorTree(const Tree& tree)
		: TreeAttribute<Tensor<T>>(tree) {
	}

	void normalize(double eps = 1e-10);

	TensorTree(const Tree& tree, function<Tensor<T>(const TensorShape&)> generator);

	void print() const;
};

typedef TensorTree<complex<double>> TensorTreecd;

typedef TensorTree<double> TensorTreed;

#endif //TENSORTREE_H
