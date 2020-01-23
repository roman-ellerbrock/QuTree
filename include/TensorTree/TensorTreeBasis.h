//
// Created by Roman Ellerbrock on 2020-01-21.
//

#ifndef TENSORTREEBASIS_H
#define TENSORTREEBASIS_H
#include "Node.h"
#include "LinearizedLeaves.h"
#include <map>

typedef vector<reference_wrapper<Node>> LinearizedNodes;

class TensorTreeBasis {
public:
	/// Default constructor
	TensorTreeBasis() = default;

	/// Default Destructor
	~TensorTreeBasis() = default;

	/// File constructor
	explicit TensorTreeBasis(istream& is);

	explicit TensorTreeBasis(const string& filename);

	/// Create Balanced Tree
	TensorTreeBasis(size_t order, size_t dim_leaves, size_t dim_nodes);

	/// Copy constructor
	TensorTreeBasis(const TensorTreeBasis& T);

	/// Move constructor
	TensorTreeBasis(TensorTreeBasis&& T) noexcept;

	/// Copy assignment operator
	TensorTreeBasis& operator=(const TensorTreeBasis& T);

	/// Move assignment operator
	TensorTreeBasis& operator=(TensorTreeBasis&& T) noexcept;

	/// Read Basis from ASCII-file
	void Read(istream& is);
	void Read(const string& filename);

	/// Write the basis to a stream in ASCII output
	void Write(ostream& os = cout) const;

	/// Re-initialize the basis from Tree
	void Update();

	/// Print out basis information
	void info(ostream& os = cout) const;

	/// number of Nodes
	size_t nTotalNodes() const { return tree.nTotalNodes(); }

	/// number of logical nodes
	size_t nNodes() const { return tree.nNodes(); }

	/// number of physical nodes
	size_t nLeaves() const { return tree.nLeaves(); }

	/// Number of states
	size_t nStates() const { return TopNode().TDim().getntensor(); }

	/// Return the reference to the next node.
	/// This routine is only used for initialization once.
	/// Please use the iterator, or MCTDHNode(i) functions to address nodes!
	AbstractNode& nextNode() { return *tree.nextNode(); }

	/// get reference to Physical Coordinate i/nPhysNodes
	Leaf& GetLeaf(size_t i);
	const Leaf& GetLeaf(size_t i) const;

	/// get reference to mctdh-node i/nmctdhNodes
	Node& GetNode(size_t i);
	const Node& GetNode(size_t i) const;

	/// get reference to the mctdh topnode
	Node& TopNode() { return linearizedNodes_.back(); }

	// @TODO: Use vector-member back()
	const Node& TopNode() const { return linearizedNodes_.back(); }

	void ReindexLeafModes(map<size_t, size_t> Map);

	/// Expand a node in the Basis
	void ExpandNode(Node& node);

	void ReplaceNode(Node& old_node, Node& new_node);

	/// Bottom-up iterator over all nodes in the mctdh-tree
	/// For top-up iteration examples refer to e.g. the density-matrix class.
	vector<reference_wrapper<Node>>::const_iterator begin() const {
		return linearizedNodes_.begin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	/// For top-up iteration examples refer to e.g. the density-matrix class.
	vector<reference_wrapper<Node>>::const_iterator end() const {
		return linearizedNodes_.end();
	}

protected:
	void LinearizeNodes();

	void LinearizeLeaves();

	/// Reference block to physical coordinates
	LinearizedLeaves linearizedLeaves_;

	/// Reference block to mctdh-nodes
	LinearizedNodes linearizedNodes_;

	/// MCTDH tree holds memory
	Node tree;
};

typedef TensorTreeBasis TTBasis;

ostream& operator<<(ostream& os, TTBasis& basis);
istream& operator<<(istream& is, TTBasis& basis);

#endif //TENSORTREEBASIS_H
