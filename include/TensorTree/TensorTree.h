//
// Created by Roman Ellerbrock on 2020-01-23.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "TreeStructuredObject.h"
#include "TensorTreeBasis.h"
#include "Tensor_Implementation.h"

template<typename T>
class TensorTree: public TreeStructuredObject<Tensor<T>> {
public:
	using TreeStructuredObject<Tensor<T>>::attributes;

	/// Default constructor without memory allocation
	TensorTree() = default;
	/// Constructor with allocation of memory
	explicit TensorTree(const TTBasis& basis);

	explicit TensorTree(istream& is);

	explicit TensorTree(const string& filename);

	/// Default destructor
	~TensorTree() = default;

	/// Create Tensors for all nodes
	virtual void Initialize(const TTBasis& basis);

	/// Generate TTs
	void Generate(mt19937& gen, const TTBasis& basis,
		bool delta_lowest = true);

	/// (File) I/O
	void Read(istream& is);

	void Write(ostream& os) const;

	void Write(const string& filename) const;

	void print(const TTBasis& basis, ostream& os = cout) const;

	/// Setters & Getters
	size_t size() const { return attributes.size(); }

protected:
	void FillBottom(Tensor<T>& Phi, const Node& node);
	void FillUpper(Tensor<T>& Phi, mt19937& gen,
		const Node& node, bool delta_lowest = true);
};

typedef TensorTree<complex<double>> TensorTreecd;

typedef TensorTree<double> TensorTreed;

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t);
template<typename T>
istream& operator>>(istream& is, TensorTree<T>& t);


#endif //TENSORTREE_H
