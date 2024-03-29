#pragma once
#include "Core/stdafx.h"
#include "AbstractNode.h"
#include "Leaf.h"
#include "NodePosition.h"


class Node
	: public AbstractNode
/**
 * \class Node
 * \ingroup TTBasis
 * \brief This class manages a node in the tree-structured TTBasis representation.
 *
 * The class holds the tensor-dimensions (TensorDim) at the current node and holds the
 * connectivity to the parent & all child nodes (parent and child(s)). The node also
 * knows wether it is a toplayer or bottomlayer-node.
 * If the Node is a bottomlayer node, it holds a Leaf. Accessing Leaf or sublying
 * Nodes (non-leaves) for incorrect nodes will stop the program.
 * Please use the TTBasis iterator to swipe over all nodes in a TTBasis. The parent and child
 * getters should only be used to access local connectivity.
 * */
{
public:
	Node();
	Node(istream& file, Node *up, const NodePosition& position);
	Node(const Node& node);
	Node(Node&& node) noexcept;
	Node& operator=(const Node& old);
	Node& operator=(Node&& old) noexcept;
	Node(const Leaf& leaf, size_t ntensor);
	~Node() override = default;

	// Initialize node
	void initialize(istream& file, Node *up, const NodePosition& position);

	// print out information to this node
	void info(ostream& os = cout) const override;
	// Write the node information
	void write(ostream& os) const override;

	// Setter & Getter
	// number of sublaying nodes
	size_t nTotalNodes() const override { return nTotalNodes_; }

	// number of sublaying mctdhNodes
	size_t nNodes() const override { return nNodes_; }

	// Number of Physical Nodes at bottom under this node
	size_t nLeaves() const override { return nLeaves_; }

	// Type of node
	int type() const override { return nodeType_; }

	// True if this node is a bottomLayer_-node, otherwise false
	bool isBottomlayer() const { return bottomLayer_; }

	// True if this node is a toplayer-node, otherwise false
	bool isToplayer() const;
	// Getter for the Physical Coordinate (only if this is a bottomLayer_ node)
	Leaf& getLeaf();
	const Leaf& getLeaf() const;
	// Getter for Parent-AbstractNode
	const Node& parent() const;
	Node& parent();
	// Return reference to the i-th child
	const Node& child(size_t i) const;
	Node& child(size_t i);

	// Returns the index in the vector of children of the parent
	// (e.g. this is the 2nd child: this getter returns "1")
	int childIdx() const { return position_.childIdx(); }
	size_t parentIdx() const { return nChildren(); }

	// Getter for the number of children of this node
	int nChildren() const { return child_.size(); }

	// Getter for the TensorDim
	TensorShape& shape() { return tensorDim_; }

	const TensorShape& shape() const { return tensorDim_; }

	// Expand one of the children node in the multilayer representation
	void expandChild(size_t i);

	// Get position_ index
	NodePosition position() const { return position_; }

	void push_back(const Node& node);

	vector<unique_ptr<AbstractNode>>::const_iterator begin() const {
		return child_.begin();
	}

	vector<unique_ptr<AbstractNode>>::const_iterator end() const {
		return child_.end();
	}

	// Danger-zone (take care what you do here!)
	// Do not access the functions below, unless you really know what you are doing!
	// Setter for the address_ of this node
	void setAddress(int newaddress) { address_ = newaddress; }

	// Getter for the address_ of this node
	int address() const { return address_; }

	// pointer to the next node in sweep
	AbstractNode *nextNode() override;

    // pointer to ne next node in an SCF-sense
    AbstractNode *nextSCFNode(AbstractNode* in) override;

	// sween for pointer to the next node in sweep. Same sweep like
	// in Uwe Manthe's fortran code
	AbstractNode *nextNodeManthe() override;
	// Move getter for unique_ptr to children
	unique_ptr<AbstractNode> downUnique(size_t i);

	// Set the upwards pointer
	void setParent(AbstractNode *up) override { parent_ = up; }

	// Replace a Child
	void replace(Node& new_child, size_t idx);
	// Update counters, position indices
	void update(const NodePosition& p) override;
	// Update position_ index
	void updatePosition(const NodePosition& p);
	// Update nTotalNodes_, nNodes_ and nLeaves_
	void updatennodes();
	// Update the TensorDim
	void updateTDim();
	// Reset all counters for the swipe
	void resetCounters();

	// Get a reference to the top-node
	Node& topNode();

protected:
	// number of nodes under the current node plus this node (n_children+1)
	size_t nTotalNodes_;
	size_t nNodes_;
	size_t nLeaves_;

	// TensorDim that belongs to this node
	TensorShape tensorDim_;

	// pointer to the upwards node
	AbstractNode *parent_;
	// vector of references to the children nodes
	vector<unique_ptr<AbstractNode>> child_;

	// reference to the last_ node that was pointed at at a sweep through the layer_
	int nextNodeNum_{0};
	size_t nextNodeNumFortran_{0};

	// position object
	NodePosition position_;
	int address_;

	// type of this node (PhysicalNode=0, GetNode=1, oSQRNode=2)
	int nodeType_;
	bool bottomLayer_;
};
