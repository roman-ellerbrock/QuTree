#pragma once
#include "stdafx.h"
#include "AbstractNode.h"
#include "Leaf.h"
#include "NodePosition.h"

/**
 * \class mctdhNode
 * \ingroup Basis
 * \brief This class manages a node in the tree-structured MCTDH-basis representation.
 *
 * The class holds the tensor-dimensions (TensorDim) at the current node and holds the
 * connectivity to neighboring parent & child nodes (Up and Down(s)). The node also 
 * knows wether it is a toplayer or bottomlayer-node.
 * If the mctdhNode is a bottomlayer node (often called a ``leave''), if holds a
 * PhysicalCoordinate. Accessing PhysCoord (leaves) or sublying mctdhNodes (non-leaves)
 * for incorrect nodes will stop the program.
 *
 * */

class Node
	: public AbstractNode {
public:
	Node();
	Node(istream& file, Node *up, const NodePosition& position);
	Node(const Node& node);
	Node(Node&& node) noexcept;
	Node& operator=(const Node& old);
	Node& operator=(Node&& old) noexcept;
	Node(const Leaf& phys, size_t ntensor);
	~Node() override = default;

	// Initialize node
	void Initialize(istream& file, Node *up, const NodePosition& position);

	// print out information to this node
	void info(ostream& os = cout) const override;
	// Write the node information
	void Write(ostream& os) const override;

	// Setter & Getter
	// number of sublaying nodes
	size_t nTotalNodes() const override { return nTotalNodes_; }

	// number of sublaying mctdhNodes
	size_t nNodes() const override { return nNodes_; }

	// Number of Physical Nodes at bottom under this node
	size_t nLeaves() const override { return nLeaves_; }

	// Type of node
	int NodeType() const override { return nodeType_; }

	// True if this node is a bottomLayer_-node, otherwise false
	bool IsBottomlayer() const { return bottomLayer_; }

	// True if this node is a toplayer-node, otherwise false
	bool IsToplayer() const;
	// Getter for the Physical Coordinate (only if this is a bottomLayer_ node)
	Leaf& PhysCoord();
	const Leaf& PhysCoord() const;
	// Getter for Parent-AbstractNode
	const Node& Up() const;
	Node& Up();
	// Return reference to the i-th child
	const Node& Down(size_t i) const;
	Node& Down(size_t i);

	// Returns the index in the vector of children of the parent
	// (e.g. this is the 2nd child: this getter returns "1")
	int ChildIdx() const { return position_.ChildIdx(); }

	// Getter for the number of children of this node
	int nChildren() const { return down_.size(); }

	// Getter for the TensorDim
	TensorDim& TDim() { return tensorDim_; }

	const TensorDim& TDim() const { return tensorDim_; }

	// Expand one of the children node in the multilayer representation
	void ExpandChild(size_t i);

	// Get position_ index
	NodePosition Position() const { return position_; }

	void push_back(const Node& node);

	vector<unique_ptr<AbstractNode>>::const_iterator begin() const {
		return down_.begin();
	}

	vector<unique_ptr<AbstractNode>>::const_iterator end() const {
		return down_.end();
	}

	// Danger-zone (take care what you do here!)
	// Do not access the functions below, unless you really know what you are doing!
	// Setter for the address_ of this node
	void SetAddress(int newaddress) { address_ = newaddress; }

	// Getter for the address_ of this node
	int Address() const { return address_; }

	// pointer to the next node in sweep
	AbstractNode *nextNode() override;
	// sween for pointer to the next node in sweep. Same sweep like
	// in Uwe Manthe's fortran code
	AbstractNode *nextNodeManthe() override;
	// Move getter for unique_ptr to children
	unique_ptr<AbstractNode> DownUnique(size_t i);
	// Set the upwards pointer
	void SetUp(AbstractNode* up) override { up_ = up; }

	// Replace a Child
	void Replace(Node& new_child, size_t idx);
	// Update counters, position indices
	void Update(const NodePosition& p) override;
	// Update position_ index
	void UpdatePosition(const NodePosition& p);
	// Update nTotalNodes_, nNodes_ and nLeaves_
	void Updatennodes();
	// Update the TensorDim
	void UpdateTDim();
	// Reset all counters for the swipe
	void ResetCounters();

	// Get a reference to the top-node
	Node& TopNode();

protected:
	// number of nodes under the current node plus this node (n_children+1)
	size_t nTotalNodes_;
	size_t nNodes_;
	size_t nLeaves_;

	// TensorDim that belongs to this node
	TensorDim tensorDim_;

	// pointer to the upwards node
//	Node *up_;
	AbstractNode *up_;
	// vector of references to the children nodes
	vector<unique_ptr<AbstractNode>> down_;

	// reference to the last node that was pointed at at a sweep through the layer_
	int nextNodeNum_;
	size_t nextNodeNumFortran_;

	// Position object
	NodePosition position_;
	int address_;

	// type of this node (PhysicalNode=0, GetNode=1, oSQRNode=2)
	int nodeType_;
	bool bottomLayer_;
};
