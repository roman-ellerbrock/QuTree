//
// Created by Roman Ellerbrock on 12/3/21.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "SubTreeAttribute.h"

template<typename T, template <typename> class Memory = polymorphic::hostMemory>
class TensorTree: public SubTreeAttribute<Tensor<T, Memory>> {
public:
	TensorTree() = default;
	~TensorTree() = default;

	explicit TensorTree(const Tree& tree,
		function<Tensor<T, Memory>(const TensorShape&)> generator = randomGen);

	explicit TensorTree(const SubTree& subTree) : SubTreeAttribute<Tensor<T,Memory>>(subTree) {
	}

	void normalize(double eps = 1e-10);

	void print() const override;
};

typedef TensorTree<complex<double>> TensorTreecd;

typedef TensorTree<double> TensorTreed;

template<typename T, template <typename> class Mem = polymorphic::hostMemory>
Tensor<T,Mem> normalize(const Tensor<T,Mem>& phi, const Edge* edge, double eps = 1e-10);

template<typename T, template <typename> class Mem = polymorphic::hostMemory>
ostream& output(ostream& os, const TensorTree<T, Mem>& A);

template<typename T, template <typename> class Mem = polymorphic::hostMemory>
ostream& operator<<(ostream& os, const TensorTree<T, Mem>& A);

#endif //TENSORTREE_H
