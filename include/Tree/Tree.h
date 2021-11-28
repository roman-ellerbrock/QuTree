//
// Created by Roman Ellerbrock on 2020-01-21.
//
#ifndef TREE_H
#define TREE_H
#include "Node.h"
#include "Edge.h"
#include "LeafArray.h"
#include <map>

typedef vector<reference_wrapper<Node>> NodeArray;

class Tree {
	/**
	 * \class Tree
	 * \ingroup Tree
	 * \brief This class manages the tensor tree tree.
	 *
	 * Tree holds and manages the tree topology
	 * and TensorShape at every node. It provides random access
	 * to nodes and forward/backward iterators
	 *
	 * Usage:
	 * Tree tree(n_leaves dim_leaves, dim_nodes); // Create close to balanced tree
	 * for (const Node& node : tree.bottomUp()) {
	 * 		// Do something - bottom-up swipe
	 * }
	 *
	 * for (const Node& node : tree.topDown()) {
	 * 		// Do something - Top-child swipe
	 * }
	 * or
	 *
	 * for (const Leaf& leaf : tree.leaves()) {
	 * 		// Do something - for every leaf
	 * }
	 */
public:
	/// Default constructor
	Tree() = default;

	/// Default Destructor
	~Tree() = default;

	explicit Tree(istream& is) { read(is); }
	explicit Tree(const string& filename) { read(filename); }

	Tree(const Tree& T);
	Tree(Tree&& T) noexcept;
	Tree& operator=(const Tree& T);
	Tree& operator=(Tree&& T) noexcept;

	/// read from ASCII-file
	void read(istream& is);
	void read(const string& filename);

	/// Write in ASCII format
	void write(ostream& os = cout) const;

	/// Print out tree information
	void info(ostream& os = cout) const;

	/// Human readable output of the tree shape
	void print(ostream& os = cout) const;

	/// number of logical nodes
	[[nodiscard]] size_t nNodes() const { return root_.nNodes(); }

	/// number of physical nodes
	[[nodiscard]] size_t nLeaves() const { return root_.nLeaves(); }

	/// Number of states
	[[nodiscard]] size_t nStates() const { return root().shape_.lastDimension(); }

	[[nodiscard]] Leaf& getLeaf(size_t mode);
	[[nodiscard]] const Leaf& getLeaf(size_t mode) const;

	/// get reference to mctdh-node i/nmctdhNodes
	[[nodiscard]] Node& getNode(size_t address);
	[[nodiscard]] const Node& getNode(size_t address) const;

	[[nodiscard]] Node& root() { return nodeArray_.back(); }

	[[nodiscard]] const Node& root() const { return nodeArray_.back(); }

//	[[nodiscard]] auto bottomUp() const;

	[[nodiscard]] auto topDown() const { return nodeArray_; }

	[[nodiscard]] const NodeArray& nodes() const { return nodeArray_; }

	[[nodiscard]] const vector<Edge>& edges() const { return edges_; }

	/// Set the root of the tree and update the Tree
	void setRoot(Node& root) {
		root_ = root;
		update();
	}

	/// Bottom-up iterator over all nodes in the tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_iterator begin() const {
		return nodeArray_.begin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_iterator end() const {
		return nodeArray_.end();
	}

	/// Top-down iterator over all nodes in the mctdh-tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_reverse_iterator rbegin() const {
		return nodeArray_.rbegin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_reverse_iterator rend() const {
		return nodeArray_.rend();
	}

protected:
	/// reinitialize from root node
	void update();

	void linearizeNodes();

	/// Reference block to leaves
	LeafArray leafArray_;

	/// Reference block to nodes
	NodeArray nodeArray_;

	vector<Edge> edges_;

	/// MCTDH tree holds memory
	Node root_;
};

ostream& operator<<(ostream& os, const Tree& tree);
istream& operator>>(istream& is, Tree& tree);

/// Assign new indices to leaves
//	void reindexLeafModes(map<size_t, size_t> Map);
/// Reset indices of leaves
//	void resetLeafModes();

/// replace a node in the tree with a new node
//	void replaceNode(Node& old_node, Node& new_node);
//	void expandNode(Node& node);


#endif //TREE_H
