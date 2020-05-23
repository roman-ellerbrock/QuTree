//
// Created by Roman Ellerbrock on 2019-07-14.
//
#ifndef MATRIXLISTTREE_H
#define MATRIXLISTTREE_H
#include "TreeOperators/TreeStructured/TreeSOP.h"
#include "TreeClasses/TensorTree.h"

typedef vector<Matrixcd> MatrixList;

class MatrixListTree: public NodeAttribute<MatrixList> {
public:
	/*!
	 * \brief This class manages multilayer H-matrices of multilayer SOP operators.
	 *
	 * \ingroup mlMatrices
	 *
	 *
	 *
	 */

	/// Default constructor
	MatrixListTree() = default;

	/// Default Destructor
	~MatrixListTree() = default;

	/// Constructor that initializes the mlHMatrices
	MatrixListTree(const TreeSOP& op, const Tree& tree)
		{ Initialize(op, tree); }

	/// Initialize the mlHMatrices
	void Initialize(const TreeSOP& op, const Tree& tree);

	/// Calculating mlHMatrices for a given wavefunction and a mNodeSOP operator
	void Calculate(const TensorTreecd& Psi, const TreeSOP& op,
		const Tree& tree);

	/// Apply a NodeSOP operator
	Tensorcd Apply(const Tensorcd& Phi, const NodeSOP& S,
		const TreeSOP& H,
		const Node& node) const;

	/// Hole-Apply a NodeSOP operator
	Tensorcd ApplyHole(const Tensorcd& Phi, const NodeSOP& S,
		const Node& node, const NodeOperator& hole_S) const;

	/// Print the ml-H-matrices
	void print(const Tree& tree, ostream& os = cout) const;

private:
	/// Calculating mlHMatrices for a TreeSOP
	MatrixList Calculate(const Tensorcd& Phi, const TreeSOP& H,
		const Node& node);

	/// Calculating mlHMatrices for a NodeSOP
	Matrixcd Calculate(const Tensorcd& Phi, const NodeSOP& S,
		const TreeSOP& H, const Node& node);

	/// Apply a lMPO
	Tensorcd ApplyNodeProductOperator(const Tensorcd& Phi,
		const NodeProductOperator& M, const TreeSOP& H,
		const Node& node) const;

	/// Apply a lMPO at an upper layer
	Tensorcd ApplyNodeProductOperator_upper(Tensorcd Phi,
		const NodeProductOperator& M, const Node& node) const;

	/// Apply a lMPO at the bottomlayer
	Tensorcd ApplyNodeProductOperator_bottom(Tensorcd Phi,
		const NodeProductOperator & M,
		const LeafOperatorLib& lib, const Node& node) const;

	/// Hole-Apply a lMPO
	Tensorcd ApplyHole(Tensorcd MPhi, const NodeProductOperator & M,
		const Node& node, const NodeOperator& hole_S) const;
};


#endif //MATRIXLISTTREE_H
