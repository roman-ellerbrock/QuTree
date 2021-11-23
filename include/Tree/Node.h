#pragma once
#include "stdafx.h"
#include "Leaf.h"
#include "NodePosition.h"
#include "Tensor/Tensor"

class Node
/**
 * \class Node
 * \ingroup Tree
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
	Node(const Node& node);
	Node(Node&& node) noexcept;
	Node& operator=(const Node& old);
	Node& operator=(Node&& old) noexcept;
	~Node() = default;

	Node(istream& file, Node *up, const NodePosition& position);
	Node(const Leaf& leaf, size_t ntensor);

	void initialize(istream& file, Node *up, const NodePosition& position);

	// print out information to this node
	void info(ostream& os = cout) const;
	// Write the node information
	void write(ostream& os) const;

	[[nodiscard]] size_t nNodes() const {
		size_t n = 1;
		for (const Node& child : child_) {
			n += child.nNodes();
		}
		return n;
	}

	[[nodiscard]] size_t nLeaves() const {
		size_t n = 0;
		for (const Node& child : child_) {
			n += child.nLeaves();
		}
		for (const Leaf& leaf : leaves_) {
			n += 1;
		}
		return n;
	}

	// True if this node is a bottomLayer_-node, otherwise false
	[[nodiscard]] bool isBottomlayer() const { return (!leaves_.empty()); }

	// True if this node is a toplayer-node, otherwise false
	[[nodiscard]] bool isToplayer() const { return (parent_ == nullptr); }

	// Getter for Parent-AbstractNode
	const Node& parent() const;
	Node& parent();

	// Returns the index in the vector of children of the parent
	// (e.g. this is the 2nd child: this getter returns "1")
	int childIdx() const { return position_.childIdx(); }
	size_t parentIdx() const { return nChildren(); }

	// Getter for the number of children of this node
	int nChildren() const { return child_.size(); }

	// Expand one of the children node in the multilayer representation
//	void expandChild(size_t i);

	// Get position_ index
	NodePosition position() const { return position_; }

	// Danger-zone (take care what you do here!)
	// Do not access the functions below, unless you really know what you are doing!
	// Setter for the address_ of this node
	void setAddress(int newaddress) { address_ = newaddress; }

	// Getter for the address_ of this node
	int address() const { return address_; }

	// pointer to the next node in sweep
	Node *nextNode();

	// sween for pointer to the next node in sweep. Same sweep like
	// in Uwe Manthe's fortran code
	Node *nextNodeManthe();

	// Update counters, position indices
	void update(const NodePosition& p);

	// Update position_ index
	void updatePosition(const NodePosition& p);

	// Reset all counters for the swipe
	void resetCounters();

	void reconnect();

	TensorShape shape_;
	int address_;

	Node *parent_;
	vector<Node> child_;
	vector<Leaf> leaves_;

	NodePosition position_;
protected:
	// reference to the last_ node that was pointed at at a sweep through the layer_
	int nextNodeNum_;
	size_t nextNodeNumFortran_;
};
