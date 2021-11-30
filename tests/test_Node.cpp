//
// Created by Roman Ellerbrock on 11/19/21.
//
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tree/Leaf.h"
#include "Tree/Node.h"
#include "Tree/LeafArray.h"

SUITE (Nodes) {

	class lf {
	public:
		lf() {
			par_ = BasisParameters({1., 0., 0., 1.,
									10, 2, 0});
			leaf_.initialize(par_);
		}

		BasisParameters par_;
		Leaf leaf_;
	};

	class small_tree: public lf {
	public:
		small_tree() {
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

	TEST_FIXTURE (lf, Leaf) {
			CHECK_EQUAL(true, leaf_.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == leaf_.api_.basis()->par_);
	}

	TEST_FIXTURE (lf, leaf_copy) {
		Leaf l2(leaf_);
			CHECK_EQUAL(true, l2.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == l2.api_.basis()->par_);
	}

	TEST_FIXTURE (lf, leaf_equal) {
		Leaf l2 = leaf_;
			CHECK_EQUAL(true, l2.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == l2.api_.basis()->par_);

		BasisParameters par2 = par_;
		par2.mode_ = par_.mode_ + 1;
			CHECK_EQUAL(false, par_ == par2);
		l2.par() = par2;
			CHECK_EQUAL(false, leaf_ == l2);
	}

	// ==================================================
	// ==== Node ========================================
	// ==================================================

	TEST (Node) {
		Node node;
			CHECK_EQUAL(1, node.nNodes());
			CHECK_EQUAL(0, node.nLeaves());
			CHECK_EQUAL(false, node.isBottomlayer());
			CHECK_EQUAL(true, node.isToplayer());
			CHECK_EQUAL(0, node.parentIdx());
	}

	TEST_FIXTURE (lf, Node_leaf) {
		Node node;
		Leaf leaf;
		leaf.initialize(par_);
		node.push_back(leaf);
			CHECK_EQUAL(1, node.nLeaves());
			CHECK_EQUAL(&node, node.leaves_[0].parent_);
			CHECK_EQUAL(true, node.leaves_[0].api_.basis()->par_ == par_);
	}

	TEST (Node_push_back) {
		Node root;
		root.push_back(Node());
			CHECK_EQUAL(2, root.nNodes());
			CHECK_EQUAL(0, root.nLeaves());
			CHECK_EQUAL(false, root.isBottomlayer());
			CHECK_EQUAL(true, root.isToplayer());
			CHECK_EQUAL(1, root.parentIdx());

		const Node& child = root.child_[0];
			CHECK_EQUAL(1, child.nNodes());
			CHECK_EQUAL(0, child.nLeaves());
			CHECK_EQUAL(false, child.isBottomlayer());
			CHECK_EQUAL(false, child.isToplayer());
			CHECK_EQUAL(0, child.parentIdx());
	}

	TEST_FIXTURE (small_tree, Node_push_back2) {
		/// singlelayer with two nodes
		/// check first child
			CHECK_EQUAL(3, root_.child_[0].nNodes());
			CHECK_EQUAL(2, root_.child_[0].nLeaves());
			CHECK_EQUAL(2, root_.child_[0].parentIdx());
		/// check second child
			CHECK_EQUAL(3, root_.child_[1].nNodes());
			CHECK_EQUAL(2, root_.child_[1].nLeaves());
			CHECK_EQUAL(2, root_.child_[1].parentIdx());
		/// check root_
			CHECK_EQUAL(7, root_.nNodes());
			CHECK_EQUAL(4, root_.nLeaves());
			CHECK_EQUAL(2, root_.parentIdx());
	}

	TEST_FIXTURE (small_tree, Node_parents) {
		/// check parent_
			CHECK_EQUAL(true, (nullptr == root_.parent_));
		const Node& inter0 = root_.child_[0];
		const Node& inter1 = root_.child_[1];
			CHECK_EQUAL(&root_, inter0.parent_);
			CHECK_EQUAL(&root_, inter1.parent_);
		const Node& bottom0 = inter0.child_[0];
		const Node& bottom1 = inter0.child_[1];
			CHECK_EQUAL(&inter0, bottom0.parent_);
			CHECK_EQUAL(&inter0, bottom1.parent_);
		const Node& bottom2 = inter1.child_[0];
		const Node& bottom3 = inter1.child_[1];
			CHECK_EQUAL(&inter1, bottom2.parent_);
			CHECK_EQUAL(&inter1, bottom3.parent_);
	}

	bool operator==(const vector<size_t>& a, const vector<size_t>& b) {
		if (a.size() != b.size()) { return false; }
		for (size_t i = 0; i < a.size(); ++i) {
			if (a[i] != b[i]) { return false; }
		}
		return true;
	}

	ostream& operator<<(ostream& os, const vector<size_t>& a) {
		for (auto x : a) {
			cout << x << " ";
		}
		return os;
	}

	ostream& operator<<(ostream& os, const vector<int>& a) {
		for (auto x : a) {
			cout << x << " ";
		}
		return os;
	}

	TEST_FIXTURE (small_tree, Node_position) {
		NodePosition proot = {};
			CHECK_EQUAL(proot, root_.position_);
		NodePosition pinter0({0});
		NodePosition pinter1({1});
		const Node& inter0 = root_.child_[0];
		const Node& inter1 = root_.child_[1];
			CHECK_EQUAL(pinter0, inter0.position_);
			CHECK_EQUAL(pinter1, inter1.position_);
		NodePosition pbottom0({0, 0});
		NodePosition pbottom1({0, 1});
		NodePosition pbottom2({1, 0});
		NodePosition pbottom3({1, 1});
		const Node& bottom0 = inter0.child_[0];
		const Node& bottom1 = inter0.child_[1];
		const Node& bottom2 = inter1.child_[0];
		const Node& bottom3 = inter1.child_[1];
			CHECK_EQUAL(pbottom0, bottom0.position_);
			CHECK_EQUAL(pbottom1, bottom1.position_);
			CHECK_EQUAL(pbottom2, bottom2.position_);
			CHECK_EQUAL(pbottom3, bottom3.position_);
	}

	TEST_FIXTURE (small_tree, Node_get) {
		{
			NodePosition pos({0, 1});
			const Node *p = root_.get(pos);
				CHECK_EQUAL(pos, p->position_);
		}

		{
			NodePosition posr({});
			const Node *r = root_.get(posr);
				CHECK_EQUAL(posr, r->position_);
		}
	}

	TEST_FIXTURE (small_tree, Node_sweep0) {
		Node *c = sweep(&root_);
			CHECK_EQUAL(&root_.child_[0], c);
		Node *b = sweep(c);
			CHECK_EQUAL(&c->child_[0], b);
		Node *b1 = sweep(b);
			CHECK_EQUAL(&c->child_[1], b1);
		Node *c1 = sweep(b1);
			CHECK_EQUAL(&root_.child_[1], c1);
		Node *b2 = sweep(c1);
			CHECK_EQUAL(&c1->child_[0], b2);
		Node *b3 = sweep(b2);
			CHECK_EQUAL(&c1->child_[1], b3);
		Node *null = sweep(b3);
			CHECK_EQUAL(true, nullptr == null);
	}

	TEST_FIXTURE (small_tree, Node_sweep) {
		size_t n = 0;
		Node *p = &root_;
		while ((p = sweep(p))) {
			n++;
		}
			CHECK_EQUAL(6, n);
	}

	TEST (Node_read) {
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
		CHECK_EQUAL(7, root.nNodes());
		CHECK_EQUAL(4, root.nLeaves());
	}

	TEST(Node_moveassign){
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
			CHECK_EQUAL(7, root.nNodes());
			CHECK_EQUAL(4, root.nLeaves());

	}

	TEST_FIXTURE (small_tree, node_copy) {
		Node node;
		node = root_;
			CHECK_EQUAL(root_.nNodes(), node.nNodes());
	}

	TEST_FIXTURE (small_tree, node_move) {
		Node tmp = root_;
		Node node = move(tmp);
			CHECK_EQUAL(root_.nNodes(), node.nNodes());
	}

	TEST_FIXTURE (small_tree, node_moveassign) {
		Node tmp = root_;
		Node node(move(tmp));
			CHECK_EQUAL(root_.nNodes(), node.nNodes());
	}

	TEST_FIXTURE (small_tree, node_copyconstr) {
		Node node(root_);
			CHECK_EQUAL(root_.nNodes(), node.nNodes());
	}

	TEST_FIXTURE (small_tree, node_equal) {
		Node root = root_;
			CHECK_EQUAL(true, root == root_);
		Node root2;
		root2.push_back(Node());
			CHECK_EQUAL(false, root2 == root_);
	}

	TEST_FIXTURE (small_tree, node_notequal) {
		Node root = root_;
			CHECK_EQUAL(false, root != root_);
		Node root2;
		root2.push_back(Node());
			CHECK_EQUAL(true, root2 != root_);
	}

	TEST (LinearizdLeaves) {
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
		LeafArray linearizedLeaves(root);
			CHECK_EQUAL(4, linearizedLeaves.size());
		for (size_t i = 0; i < 4; ++i) {
				CHECK_EQUAL(i, linearizedLeaves[i].api_.basis()->par_.mode_);
		}
	}
}
