//
// Created by Roman Ellerbrock on 2019-04-22.
//

#ifndef MCTDH_TREE_H
#define MCTDH_TREE_H
#include "TreeNode.h"
#include <map>

TreeNode Group(const TreeNode& node, const NodeContent& c) {
	TreeNode parent(c);
	for (size_t k = 0; k < c.f; ++k) {
		parent.push_back(node);
	}
	return parent;
}

class TensorDimTree {
public:
	TensorDimTree() = default;
	~TensorDimTree() = default;

	TensorDimTree(size_t n_layers, size_t n, size_t N, size_t d) {
		/// Create a d-ary balanced tree with n spf
		TreeNode leaf(N);
		TreeNode bottom(n);
		TreeNode node(n);
		NodeContent c(n, 2);
		for (size_t l = 0; l < n_layers; ++l) {
			node = Group(node, c);
		}
		root_ = node;
		root_.content_.n = 1;
		root_.MakeRoot();
		size_t mode = 0;
		root_.SetPhysicalNodes(mode);

		Update();
	}

	void Update() {
		TreeNode *next = root_.NextNode();
		nodes_.clear();
		size_t idx = 0;
		while (next != &root_) {
			if (next->IsLeaf()) {
				phys_nodes_.emplace_back(reference_wrapper<TreeNode>(*next));
			} else {
				nodes_.emplace_back(reference_wrapper<TreeNode>(*next));
			}
			next->SetAddress(idx);
			idx++;
			next = root_.NextNode();
		}
		nodes_.emplace_back(reference_wrapper<TreeNode>(*next));

		root_.UpdateTensorDim();
	}

	void print(ostream& os = cout)const {
		for (const TreeNode& node : *this) {
			::print(node.GetPath(), os);
			node.TDim().print(os);
			cout << endl;
		}
	}

	vector<reference_wrapper<TreeNode>>::const_iterator begin() const{ return nodes_.begin(); }
	vector<reference_wrapper<TreeNode>>::const_iterator end() const{ return nodes_.end(); }

	size_t nLeaves()const { return phys_nodes_.size(); }
	const TreeNode& Leaf(size_t k)const { return phys_nodes_[k]; }
	TreeNode& Leaf(size_t k) { return phys_nodes_[k]; }
private:
	TreeNode root_;
	vector<reference_wrapper<TreeNode>> nodes_;
	vector<reference_wrapper<TreeNode>> phys_nodes_;
//	map<size_t, reference_wrapper<TreeNode>> phys_nodes_;
};

#endif //MCTDH_TREE_H
