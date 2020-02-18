//
// Created by Roman Ellerbrock on 2020-01-21.
//

#ifndef TREE_H
#define TREE_H
#include "Node.h"
#include "LinearizedLeaves.h"
#include <map>

typedef vector<reference_wrapper<Node>> LinearizedNodes;

class Tree {
	/**
	 * \class Tree
	 * \ingroup Tree
	 * \brief This class manages the tensor tree tree.
	 *
	 * TensorTreeBasis (TTBasis) holds and manages the tree structure
	 * and holds tensor dimensions at every node. It provides iterators
	 * for swiping over every node in a tree. For a bottom-up swipe though
	 * the tree use the iterator; for a top-down swipe use a regular for loop
	 * and get Nodes via GetNode(i).
	 *
	 * Usage:
	 * Tree tree(n_leaves dim_leaves, dim_nodes); // Create close to balanced tree
	 * for (const Node& node : tree) {
	 * 		// Do something - bottom-up swipe
	 * }
	 *
	 * for (int i = tree.nNodes() - 1; i > 0; --i) {
	 * 		const Node& node = tree.GetNode(i);
	 * 		// Do something - Top-Down swipe
	 * }
	 * or
	 * for (auto it = tree.rbegin(); it != tree.rend(); it++) {
	 *
	 * }
	 *
	 * for (size_t l = 0; l < nLeaves; ++l) {
	 * 		const Leaf& leaf = GetLeaf(l);
	 * 		// Do something - for every leaf
	 * }
	 */
public:
	/// Default constructor
	Tree() = default;

	/// Default Destructor
	~Tree() = default;

	/// File constructor
	explicit Tree(istream& is);

	/// Stream constructor
	explicit Tree(const string& filename);

	/// Copy constructor
	Tree(const Tree& T);

	/// Move constructor
	Tree(Tree&& T) noexcept;

	/// Copy assignment operator
	Tree& operator=(const Tree& T);

	/// Move assignment operator
	Tree& operator=(Tree&& T) noexcept;

	/// Read Basis from ASCII-file
	void Read(istream& is);
	void Read(const string& filename);

	/// Write the tree to a stream in ASCII output
	void Write(ostream& os = cout) const;

	/// Re-initialize the tree from Tree
	void Update();

	/// Print out tree information
	void info(ostream& os = cout) const;

	/// number of Nodes
	size_t nTotalNodes() const { return root_.nTotalNodes(); }

	/// number of logical nodes
	size_t nNodes() const { return root_.nNodes(); }

	/// number of physical nodes
	size_t nLeaves() const { return root_.nLeaves(); }

	/// Number of states
	size_t nStates() const { return TopNode().TDim().LastActive(); }

	/// Return the reference to the next node.
	/// This routine is only used for initialization once.
	/// Please use the iterator, or MCTDHNode(i) functions to address nodes!
	AbstractNode& nextNode() { return *root_.nextNode(); }

	/// get reference to Physical Coordinate i/nPhysNodes
	Leaf& GetLeaf(size_t i);
	const Leaf& GetLeaf(size_t i) const;

	/// get reference to mctdh-node i/nmctdhNodes
	Node& GetNode(size_t i);
	const Node& GetNode(size_t i) const;

	/// get reference to the mctdh topnode
	Node& TopNode() { return linearizedNodes_.back(); }

	const Node& TopNode() const { return linearizedNodes_.back(); }

	/// Assign new indices to leaves
	void ReindexLeafModes(map<size_t, size_t> Map);

	/// Reset indices of leaves
	void ResetLeafModes();

	/// Expand a node in the Basis
	void ExpandNode(Node& node);

	/// Replace a node in the tree with a new node
	void ReplaceNode(Node& old_node, Node& new_node);

	/// Set the root of the tree and update the TreeShape
	void SetRoot(Node& root) { root_ = root;
		root_.UpdatePosition(NodePosition());
		Update();
	}

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

	/// Top-down iterator over all nodes in the mctdh-tree
	vector<reference_wrapper<Node>>::const_reverse_iterator  rbegin() const {
		return linearizedNodes_.rbegin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	vector<reference_wrapper<Node>>::const_reverse_iterator  rend() const {
		return linearizedNodes_.rend();
	}
	/// Check whether TensorTreeBasis is working correctly
	bool IsWorking();

	/// Human readable output of the tree shape
	void print() const;

protected:
	void LinearizeNodes();

	void LinearizeLeaves();

	/// Reference block to physical coordinates
	LinearizedLeaves linearizedLeaves_;

	/// Reference block to mctdh-nodes
	LinearizedNodes linearizedNodes_;

	/// MCTDH tree holds memory
	Node root_;
};

ostream& operator<<(ostream& os, const Tree& tree);
istream& operator>>(istream& is, Tree& tree);

#endif //TREE_H
