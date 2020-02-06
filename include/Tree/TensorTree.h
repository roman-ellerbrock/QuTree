//
// Created by Roman Ellerbrock on 2020-01-23.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "Tree/TreeStructuredObject.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "Core/Tensor_Implementation.h"

template<typename T>
class TensorTree: public TreeStructuredObject<Tensor<T>>
	/**
	 * \class TensorTree
	 * \ingroup Tree
	 * \brief This class represents tensor trees.
	 *
	 * Usage:
	 * TRBasis basis(12, 2, 2)
	 * // Create Tensor with Zero-entry tensors at every node
	 * TensorTreecd Psi(basis);
	 * Psi.Write("filename.dat");
	 * TensorTreecd Chi("filename.dat");
	 */
{
public:
	using TreeStructuredObject<Tensor<T>>::attributes_;

	/// Default constructor without memory allocation
	TensorTree() = default;
	/// Constructor with allocation of memory
	explicit TensorTree(const TTBasis& basis);

	/// Construct TensorTree from stream
	explicit TensorTree(istream& is);

	/// Construct TensorTree from file
	explicit TensorTree(const string& filename);

	/// Create tensor tree and occupy the coefficients
	TensorTree(const TTBasis& basis,
		mt19937& gen, bool delta_lowest = true);

	/// Default destructor
	~TensorTree() = default;

	/// Create Tensors for all nodes
	virtual void Initialize(const TTBasis& basis);

	/// Generate TTs
	void Generate(const TTBasis& basis,
		mt19937& gen, bool delta_lowest = true);

	/// (File) I/O
	/// Read TensorTree from stream (binary format)
	void Read(istream& is);

	/// Write TensorTree to stream (binary format)
	void Write(ostream& os) const;

	/// Write TensorTree to file (binary format)
	void Write(const string& filename) const;

	/// Print info in human readable format
	void print(const TTBasis& basis, ostream& os = cout) const;

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
