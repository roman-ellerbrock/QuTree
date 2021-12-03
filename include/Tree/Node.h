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

	Node(istream& file, Node *up,
		const NodePosition& position);

	void write(ostream& os) const;

	void info(ostream& os = cout) const;

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
		n += leaves_.size();
		return n;
	}

	[[nodiscard]] bool isBottomlayer() const { return (!leaves_.empty()); }

	[[nodiscard]] bool isToplayer() const { return (parent_ == nullptr); }

	void push_back(const Node& child) {
		child_.push_back(child);
		reconnect();
	}

	void push_back(const Leaf& leaf) {
		leaves_.push_back(leaf);
		reconnect();
	}

	[[nodiscard]] const Node& parent() const;
	Node& parent();

	// Returns the index in the vector of children of the parent
	// (e.g. this is the 2nd child: this getter returns "1")
	[[nodiscard]] size_t childIdx() const { return position_.childIdx(); }

	[[nodiscard]] size_t parentIdx() const {
		return nChildren() + leaves_.size();
	}

	[[nodiscard]] size_t nChildren() const { return child_.size(); }

	[[nodiscard]] NodePosition position() const { return position_; }

	void updatePosition(const NodePosition& p);
	void reconnect();

	[[nodiscard]] const Node *get(NodePosition p) const;

	size_t idx(const Node& node) const;

	TensorShape shape_;
	int address_;

	Node *parent_;
	vector<Node> child_;
	vector<Leaf> leaves_;

	NodePosition position_;
};

Node* sweep(Node* last);
Node readNode(istream& file, Node *up = nullptr, const NodePosition& position = NodePosition());

bool operator==(const Node& a, const Node& b);
bool operator!=(const Node& a, const Node& b);

ostream& operator<<(ostream& os, const Node& node);