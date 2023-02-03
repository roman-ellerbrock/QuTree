//
// Created by Roman Ellerbrock on 2020-01-21.
//

#ifndef TREE_H
#define TREE_H
#include "Node.h"
#include "Edge.h"
#include "LinearizedLeaves.h"
#include <map>
#include <numeric>

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
	 * 		// Do something - Top-child swipe
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

	/// read Basis from ASCII-file
	void read(istream& is);
	void read(const string& filename);

	/// Write the tree to a stream in ASCII output
	void write(ostream& os = cout) const;

	/// re-initialize the tree from Tree
	void update();

	/// Print out tree information
	void info(ostream& os = cout) const;

	/// number of Nodes
	[[nodiscard]] size_t nTotalNodes() const { return root_.nTotalNodes(); }

	/// number of logical nodes
	[[nodiscard]] size_t nNodes() const { return root_.nNodes(); }

	/// number of physical nodes
	[[nodiscard]] size_t nLeaves() const { return root_.nLeaves(); }

	/// Number of states
	[[nodiscard]] size_t nStates() const { return topNode().shape().lastDimension(); }

	/// Return the reference to the next node.
	/// This routine is only used for initialization once.
	/// Please use the iterator, or mctdhNode(i) functions to address nodes!
	AbstractNode& nextNode() { return *root_.nextNode(); }

	/// get reference to Physical Coordinate i/nPhysNodes
	Leaf& getLeaf(size_t i);
	[[nodiscard]] const Leaf& getLeaf(size_t i) const;

	/// get reference to mctdh-node i/nmctdhNodes
	Node& getNode(size_t i);
	[[nodiscard]] const Node& getNode(size_t i) const;

	/// get reference to the mctdh topnode
	Node& topNode() { return linearizedNodes_.back(); }

	[[nodiscard]] const Node& topNode() const { return linearizedNodes_.back(); }

	/// Assign new indices to leaves
	void reindexLeafModes(map<size_t, size_t> Map);

	/// Reset indices of leaves
	void resetLeafModes();

	/// Expand a node in the Basis
	void expandNode(Node& node);

	/// replace a node in the tree with a new node
	void replaceNode(Node& old_node, Node& new_node);

	/// Set the root of the tree and update the TreeShape
	void setRoot(Node& root) {
		root_ = root;
		root_.updatePosition(NodePosition());
		update();
	}

	/// Bottom-up iterator over all nodes in the mctdh-tree
	/// For top-up iteration examples refer to e.g. the density-matrix class.
	[[nodiscard]] vector<reference_wrapper<Node>>::const_iterator begin() const {
		return linearizedNodes_.begin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	/// For top-up iteration examples refer to e.g. the density-matrix class.
	[[nodiscard]] vector<reference_wrapper<Node>>::const_iterator end() const {
		return linearizedNodes_.end();
	}

	/// Top-down iterator over all nodes in the mctdh-tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_reverse_iterator rbegin() const {
		return linearizedNodes_.rbegin();
	}

	/// Bottom-up const iterator over all nodes in the mctdh-tree
	[[nodiscard]] vector<reference_wrapper<Node>>::const_reverse_iterator rend() const {
		return linearizedNodes_.rend();
	}

	/// Check whether TensorTreeBasis is working correctly
	bool isWorking();

	/// Human readable output of the tree shape
	void print(ostream& os = cout) const;

	[[nodiscard]] const vector<Edge>& edges() const { return edges_; }

	[[nodiscard]] const LinearizedNodes& nodes() const { return linearizedNodes_; }

	/// return range of 1 to n_leaves
	[[nodiscard]] vector<size_t> rangeLeaves() const {
		vector<size_t> vec(nLeaves());
		std::iota(vec.begin(), vec.end(), 0);
		return vec;
	}

	[[nodiscard]] vector<const Node*> neighbors(const Node& from, int hole = -1) const;

protected:
	void linearizeNodes();

	void linearizeLeaves();

	/// Reference block to physical coordinates
	LinearizedLeaves linearizedLeaves_;

	/// Reference block to mctdh-nodes
	LinearizedNodes linearizedNodes_;

	vector<Edge> edges_;

	/// MCTDH tree holds memory
	Node root_;
};

ostream& operator<<(ostream& os, const Tree& tree);
istream& operator>>(istream& is, Tree& tree);

#endif //TREE_H
