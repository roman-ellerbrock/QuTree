//
// Created by Roman Ellerbrock on 11/19/21.
//
#include <gtest/gtest.h>
#include <iostream>
#include "Tree/Leaf.h"
#include "Tree/Node.h"
#include "Tree/LeafArray.h"

class lf : public ::testing::Test
{
public:
	lf()
	{
		par_ = BasisParameters({1., 0., 0., 1.,
								10, 2, 0});
		leaf_.initialize(par_);
	}

	BasisParameters par_;
	Leaf leaf_;
};

class small_tree : public lf
{
public:
	small_tree()
	{
		/** the small tree:
		 *	1	-2
		 *	 	1	-2
		 *	 		1	-1
		 *	 			1	0	0
		 *	 		1	-1
		 *	 			1	0	0
		 *	 	1	-2
		 *	 		1	-1
		 *	 			1	0	0
		 *	 		1	-1
		 *	 			1	0	0
		 */
		Node inter;
		Node bottom;
		Leaf leaf;
		leaf.initialize(par_);
		bottom.push_back(leaf);
		inter.push_back(bottom);
		inter.push_back(bottom);
		root_.push_back(inter);
		root_.push_back(inter);
	}

	Node root_;
};

// ==================================================
// ==== Leaf ========================================
// ==================================================

TEST_F(lf, Leaf)
{
	EXPECT_EQ(true, leaf_.parent_ == nullptr);
	EXPECT_EQ(true, par_ == leaf_.basis_.ptr()->par_);
}

TEST_F(lf, leaf_copy)
{
	Leaf l2(leaf_);
	EXPECT_EQ(true, l2.parent_ == nullptr);
	EXPECT_EQ(true, par_ == l2.basis_.ptr()->par_);
}

TEST_F(lf, leaf_equal)
{
	Leaf l2 = leaf_;
	EXPECT_EQ(true, l2.parent_ == nullptr);
	EXPECT_EQ(true, par_ == l2.basis_.ptr()->par_);

	BasisParameters par2 = par_;
	par2.mode_ = par_.mode_ + 1;
	EXPECT_EQ(false, par_ == par2);
	l2.par() = par2;
	EXPECT_EQ(false, leaf_ == l2);
}

// ==================================================
// ==== Node ========================================
// ==================================================

TEST(Node, Node)
{
	Node node;
	EXPECT_EQ(1, node.nNodes());
	EXPECT_EQ(0, node.nLeaves());
	EXPECT_EQ(false, node.isBottomlayer());
	EXPECT_EQ(true, node.isToplayer());
	EXPECT_EQ(0, node.parentIdx());
}

TEST_F(lf, Node_leaf)
{
	Node node;
	Leaf leaf;
	leaf.initialize(par_);
	node.push_back(leaf);
	EXPECT_EQ(1, node.nLeaves());
	EXPECT_EQ(&node, node.leaves_[0].parent_);
	EXPECT_EQ(true, node.leaves_[0].basis_.ptr()->par_ == par_);
}

TEST(Node, Node_push_back)
{
	Node root;
	root.push_back(Node());
	EXPECT_EQ(2, root.nNodes());
	EXPECT_EQ(0, root.nLeaves());
	EXPECT_EQ(false, root.isBottomlayer());
	EXPECT_EQ(true, root.isToplayer());
	EXPECT_EQ(1, root.parentIdx());

	const Node &child = root.child_[0];
	EXPECT_EQ(1, child.nNodes());
	EXPECT_EQ(0, child.nLeaves());
	EXPECT_EQ(false, child.isBottomlayer());
	EXPECT_EQ(false, child.isToplayer());
	EXPECT_EQ(0, child.parentIdx());
}

TEST_F(small_tree, Node_push_back2)
{
	/// singlelayer with two nodes
	/// check first child
	EXPECT_EQ(3, root_.child_[0].nNodes());
	EXPECT_EQ(2, root_.child_[0].nLeaves());
	EXPECT_EQ(2, root_.child_[0].parentIdx());
	/// check second child
	EXPECT_EQ(3, root_.child_[1].nNodes());
	EXPECT_EQ(2, root_.child_[1].nLeaves());
	EXPECT_EQ(2, root_.child_[1].parentIdx());
	/// check root_
	EXPECT_EQ(7, root_.nNodes());
	EXPECT_EQ(4, root_.nLeaves());
	EXPECT_EQ(2, root_.parentIdx());
}

TEST_F(small_tree, Node_parents)
{
	/// check parent_
	EXPECT_EQ(true, (nullptr == root_.parent_));
	const Node &inter0 = root_.child_[0];
	const Node &inter1 = root_.child_[1];
	EXPECT_EQ(&root_, inter0.parent_);
	EXPECT_EQ(&root_, inter1.parent_);
	const Node &bottom0 = inter0.child_[0];
	const Node &bottom1 = inter0.child_[1];
	EXPECT_EQ(&inter0, bottom0.parent_);
	EXPECT_EQ(&inter0, bottom1.parent_);
	const Node &bottom2 = inter1.child_[0];
	const Node &bottom3 = inter1.child_[1];
	EXPECT_EQ(&inter1, bottom2.parent_);
	EXPECT_EQ(&inter1, bottom3.parent_);
}

bool operator==(const vector<size_t> &a, const vector<size_t> &b)
{
	if (a.size() != b.size())
	{
		return false;
	}
	for (size_t i = 0; i < a.size(); ++i)
	{
		if (a[i] != b[i])
		{
			return false;
		}
	}
	return true;
}

ostream &operator<<(ostream &os, const vector<size_t> &a)
{
	for (auto x : a)
	{
		cout << x << " ";
	}
	return os;
}

ostream &operator<<(ostream &os, const vector<int> &a)
{
	for (auto x : a)
	{
		cout << x << " ";
	}
	return os;
}

TEST_F(small_tree, Node_position)
{
	NodePosition proot = {};
	EXPECT_EQ(proot, root_.position_);
	NodePosition pinter0({0});
	NodePosition pinter1({1});
	const Node &inter0 = root_.child_[0];
	const Node &inter1 = root_.child_[1];
	EXPECT_EQ(pinter0, inter0.position_);
	EXPECT_EQ(pinter1, inter1.position_);
	NodePosition pbottom0({0, 0});
	NodePosition pbottom1({0, 1});
	NodePosition pbottom2({1, 0});
	NodePosition pbottom3({1, 1});
	const Node &bottom0 = inter0.child_[0];
	const Node &bottom1 = inter0.child_[1];
	const Node &bottom2 = inter1.child_[0];
	const Node &bottom3 = inter1.child_[1];
	EXPECT_EQ(pbottom0, bottom0.position_);
	EXPECT_EQ(pbottom1, bottom1.position_);
	EXPECT_EQ(pbottom2, bottom2.position_);
	EXPECT_EQ(pbottom3, bottom3.position_);
}

TEST_F(small_tree, Node_get)
{
	{
		NodePosition pos({0, 1});
		const Node *p = root_.get(pos);
		EXPECT_EQ(pos, p->position_);
	}

	{
		NodePosition posr({});
		const Node *r = root_.get(posr);
		EXPECT_EQ(posr, r->position_);
	}
}

TEST_F(small_tree, Node_sweep0)
{
	Node *c = sweep(&root_);
	EXPECT_EQ(&root_.child_[0], c);
	Node *b = sweep(c);
	EXPECT_EQ(&c->child_[0], b);
	Node *b1 = sweep(b);
	EXPECT_EQ(&c->child_[1], b1);
	Node *c1 = sweep(b1);
	EXPECT_EQ(&root_.child_[1], c1);
	Node *b2 = sweep(c1);
	EXPECT_EQ(&c1->child_[0], b2);
	Node *b3 = sweep(b2);
	EXPECT_EQ(&c1->child_[1], b3);
	Node *null = sweep(b3);
	EXPECT_EQ(true, nullptr == null);
}

TEST_F(small_tree, Node_sweep)
{
	size_t n = 0;
	Node *p = &root_;
	while ((p = sweep(p)))
	{
		n++;
	}
	EXPECT_EQ(6, n);
}

TEST(Node, Node_read)
{
	string file("1	-2"
				"	 	1	-2"
				"	 		1	-1"
				" 				1	0	0"
				"	 		1	-1"
				"	 			1	0	1"
				"	 	1	-2"
				"	 		1	-1"
				"	 			1	0	2"
				"	 		1	-1"
				"	 			1	0	3");
	stringstream is(file);
	Node root = readNode(is);
	EXPECT_EQ(7, root.nNodes());
	EXPECT_EQ(4, root.nLeaves());
}

TEST(Node, Node_moveassign)
{
	string file("1	-2"
				"	 	1	-2"
				"	 		1	-1"
				" 				1	0	0"
				"	 		1	-1"
				"	 			1	0	1"
				"	 	1	-2"
				"	 		1	-1"
				"	 			1	0	2"
				"	 		1	-1"
				"	 			1	0	3");
	Node root;
	stringstream is(file);
	root = readNode(is);
	EXPECT_EQ(7, root.nNodes());
	EXPECT_EQ(4, root.nLeaves());
}

TEST_F(small_tree, node_copy)
{
	Node node;
	node = root_;
	EXPECT_EQ(root_.nNodes(), node.nNodes());
}

TEST_F(small_tree, node_move)
{
	Node tmp = root_;
	Node node = move(tmp);
	EXPECT_EQ(root_.nNodes(), node.nNodes());
}

TEST_F(small_tree, node_moveassign)
{
	Node tmp = root_;
	Node node(move(tmp));
	EXPECT_EQ(root_.nNodes(), node.nNodes());
}

TEST_F(small_tree, node_copyconstr)
{
	Node node(root_);
	EXPECT_EQ(root_.nNodes(), node.nNodes());
}

TEST_F(small_tree, node_equal)
{
	Node root = root_;
	EXPECT_EQ(true, root == root_);
	Node root2;
	root2.push_back(Node());
	EXPECT_EQ(false, root2 == root_);
}

TEST_F(small_tree, node_notequal)
{
	Node root = root_;
	EXPECT_EQ(false, root != root_);
	Node root2;
	root2.push_back(Node());
	EXPECT_EQ(true, root2 != root_);
}

// ==================================================
// ==== LeafArray =============================
// ==================================================

TEST(Node, LeafArray)
{
	string file("1	-2"
				"	 	1	-2"
				"	 		1	-1"
				" 				1	0	0"
				"	 		1	-1"
				"	 			1	0	1"
				"	 	1	-2"
				"	 		1	-1"
				"	 			1	0	2"
				"	 		1	-1"
				"	 			1	0	3");
	stringstream is(file);
	Node root = readNode(is);
	LeafArray leafarray(root);
	EXPECT_EQ(4, leafarray.size());
	for (size_t i = 0; i < 4; ++i)
	{
		EXPECT_EQ(i, leafarray[i].basis_.ptr()->par_.mode_);
	}
}
