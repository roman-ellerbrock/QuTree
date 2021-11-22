#pragma once
#include "stdafx.h"
#include "NodePosition.h"
#include "PrimitiveBasis/PrimitiveBasis.h"

class Node;

class Leaf
	/**
	 * \class Leaf
	 * \ingroup Tree
	 * \brief This class represents the dangling edges in a tree.
	 *
	 * A Leaf is a dangling edge of a TTN in a Tree.
	 * in a tree of a TTBasis. Leaves contain abstract class pointers
	 * to PrimitiveBasis which provides the interface to the problem
	 * under consideration.
	 */
{
public:
	Leaf();

//	Leaf(istream& file, Node *up, NodePosition position);
//	Leaf(size_t dim, size_t mode, size_t type, size_t subtype,
//		PhysPar par);
//	Leaf(const Leaf&);
//	Leaf(Leaf&&) = default;
//	Leaf& operator=(const Leaf&);
//	Leaf& operator=(Leaf&&) = default;
	~Leaf() = default;

	void CreatePrimitiveBasis(size_t type, size_t subtype, size_t dim);

	void info(ostream& os = cout) const ;
	void write(ostream& os = cout) const ;

	size_t nNodes() const { return 0; }

	size_t nLeaves() const { return 1; }

	// Getter & Setter
	size_t mode() const { return mode_; }

	int& mode() { return mode_; }

	size_t Type() const { return type_; }

	size_t dim() const { return dim_; }

	int type() const { return nodeType_; }

	LeafInterface& interface() { return *interface_; }
	const LeafInterface& interface() const { return *interface_; }

	// This is not a GetNode& to avoid circular dependencies
	Node& parent() const { return *parent_; };

	void setPar(PhysPar par) { par_ = par; }

	PhysPar par() const { return par_; }

	double omega() const { return par_.omega(); }

	double r0() const { return par_.r0(); }

	double wfr0() const { return par_.wfr0(); }

	double wfOmega() const { return par_.wfOmega(); }

	void update(const NodePosition& p);

	void updatePosition(const NodePosition& p);

	// Danger zone
	void setParent(Node *node) { parent_ = node; }

	int dim_, type_, mode_;
	int subType_;
	int nodeType_;
	// This is not a GetNode& to avoid circular dependencies
	Node *parent_;
	PhysPar par_;
	NodePosition position_;
protected:
	unique_ptr<LeafInterface> interface_;
};

