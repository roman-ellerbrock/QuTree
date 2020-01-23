//
// Created by Roman Ellerbrock on 2020-01-23.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "TreeStructuredObject.h"
#include "TensorTreeBasis.h"
#include "Tensor_Implementation.h"

template <typename T>
class TensorTree : public TreeStructuredObject<Tensor<T>> {
public:
	using TreeStructuredObject<Tensor<T>>::attributes;

	/// Default constructor without memory allocation
	TensorTree() = default;
	/// Constructor with allocation of memory
	explicit TensorTree(const TTBasis& basis);
	/// Default destructor
	~TensorTree() = default;

	virtual void Initialize(const TTBasis& basis);

	void Read(istream& is);

	void Write(ostream& os) const;

};

typedef TensorTree<complex<double>> TensorTreecd;
typedef TensorTree<double> TensorTreed;

template <typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t);
template <typename T>
istream& operator>>(istream& is, TensorTree<T>& t);

#endif //TENSORTREE_H
