//
// Created by Roman Ellerbrock on 2020-01-23.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "NodeAttribute.h"
#include "TreeShape/Tree.h"
#include "Core/Tensor_Implementation.h"

template<typename T>
class TensorTree:
	public NodeAttribute<Tensor<T>>
	/**
	 * \class TensorTree
	 * \ingroup Tree
	 * \brief This class represents tensor trees.
	 *
	 * Usage:
	 * TRBasis tree(12, 2, 2)
	 * // Create Tensor with Zero-entry tensors at every node
	 * TensorTreecd Psi(tree);
	 * Psi.Write("filename.dat");
	 * TensorTreecd Chi("filename.dat");
	 */
{
public:
	using NodeAttribute<Tensor<T>>::attributes_;

	/// Default constructor without memory allocation
	TensorTree() = default;
	/// Constructor with allocation of memory
	explicit TensorTree(const Tree& tree);

	/// Construct TensorTree from stream
	explicit TensorTree(istream& is);

	/// Construct TensorTree from file
	explicit TensorTree(const string& filename);

	/// Create tensor tree and occupy the coefficients
	TensorTree(const Tree& tree,
		std::mt19937& gen, bool delta_lowest = true);

	/// Default destructor
	~TensorTree() = default;

	/// Create Tensors for all nodes
	virtual void Initialize(const Tree& tree);

	/// Generate TTs
	void Generate(const Tree& tree,
		std::mt19937& gen, bool delta_lowest = true);

	/// (File) I/O
	/// Read TensorTree from stream (binary format)
	void Read(istream& is);

	/// Write TensorTree to stream (binary format)
	void Write(ostream& os) const;

	/// Write TensorTree to file (binary format)
	void Write(const string& filename) const;

	/// Print info in human readable format
	void print(const Tree& tree, ostream& os = cout) const;

protected:
	void FillBottom(Tensor<T>& Phi, const Node& node);
	void FillUpper(Tensor<T>& Phi, std::mt19937& gen,
		const Node& node, bool delta_lowest = true);
};

typedef TensorTree<complex<double>> TensorTreecd;

typedef TensorTree<double> TensorTreed;

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t);
template<typename T>
istream& operator>>(istream& is, TensorTree<T>& t);


#endif //TENSORTREE_H
