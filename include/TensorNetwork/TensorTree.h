//
// Created by Roman Ellerbrock on 12/3/21.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "TreeAttribute.h"

template<typename T>
class TensorTree : public TreeAttribute<Tensor<T>> {
public:
	TensorTree() = default;
	~TensorTree() = default;

	explicit TensorTree(const Tree& tree)
		: TreeAttribute<Tensor<T>>(tree) {

	}

	void normalize();

	TensorTree(const Tree& tree, function<Tensor<T>(const TensorShape&)> generator);
};

typedef TensorTree<complex<double>> TensorTreecd;
typedef TensorTree<double> TensorTreed;

template <typename T>
[[nodiscard]] TensorTree<T> sizedTensorTree(const Tree& tree);

constexpr auto sizedTensorTreecd = sizedTensorTree<complex<double>>;
constexpr auto sizedTensorTreed = sizedTensorTree<double>;

template <typename T>
TensorTree<T> randomTensorTree(const Tree& tree);

constexpr auto randomTensorTreecd = randomTensorTree<complex<double>>;
constexpr auto randomTensorTreed  = randomTensorTree<double>;



#endif //TENSORTREE_H
