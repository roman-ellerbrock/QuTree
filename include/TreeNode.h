//
// Created by Roman Ellerbrock on 2019-04-22.
//

#ifndef MCTDH_TREENODE_H
#define MCTDH_TREENODE_H
#include <utility>
#include <vector>
#include <iostream>
#include "TensorDim.h"

typedef vector<size_t> Path;

void print(const Path& p) {
	cout << "{ ";
	if (p.empty()) {
		cout << "{ }" << endl;
		return;
	}
	cout << p[0];
	for (size_t i = 1; i < p.size(); ++i) {
		cout << ", " << p[i];
	}
	cout << " }" << endl;
}

class NodeContent {
public:
	NodeContent()
		: f(0), n(0) {}

	explicit NodeContent(size_t n_)
		: f(0), n(n_) {}

	~NodeContent() = default;

	void print(size_t indent = 0, std::ostream& os = std::cout) const {
		for (size_t i = 0; i < indent; ++i) {
			cout << "\t";
		}
		os << f << " " << n << std::endl;
	}

	size_t f, n;
};

class TreeNode {
public:

	/**
	 * A simple tree-node class.
	 */

	/// Rule of five-section (constructors & destructors)
	TreeNode()
		: up(nullptr), nextnodenr_(0), address_(0) {}

	explicit TreeNode(size_t n)
		: up(nullptr), nextnodenr_(0), content_(n), address_(0) {}

	explicit TreeNode(NodeContent con)
		: up(nullptr), nextnodenr_(0), content_(con), address_(0) {}

	~TreeNode() = default;

	/// Copy constructor
	TreeNode(const TreeNode& node)
		: up(node.up), content_(node.content_),
		  nextnodenr_(node.nextnodenr_), address_(node.address_),
		  tdim_(node.tdim_) {

		for (const TreeNode *child : node.down_) {
			down_.push_back(new TreeNode(*child));
		}
		for (TreeNode *child : down_) {
			child->up = this;
		}
	}

	/// Move constructor
	TreeNode(TreeNode&& node) noexcept
		: up(node.up), down_(move(node.down_)), content_(node.content_),
		  nextnodenr_(node.nextnodenr_), address_(node.address_),
		  tdim_(node.tdim_) {

		for (TreeNode *child : down_) {
			child->up = this;
		}
	}

	/// Copy asignment operator
	TreeNode& operator=(const TreeNode& old) {
		TreeNode node(old);
		*this = std::move(node);
		return *this;
	}

	/// Move asignment operator
	TreeNode& operator=(TreeNode&& old) noexcept {
		if (this == &old) {
			return *this;
		}

		up = old.up;
		down_ = move(old.down_);
		content_ = old.content_;
		nextnodenr_ = old.nextnodenr_;
		address_ = old.address_;
		tdim_ = old.tdim_;

		for (TreeNode *child : down_) {
			child->up = this;
		}
		return *this;
	}

	/// Member functions
	TreeNode *NextNode() {
		TreeNode *result;
		if (nextnodenr_ < down_.size()) {
			result = down_[nextnodenr_]->NextNode();
			if (result == down_[nextnodenr_]) {
				++nextnodenr_;
			}
		} else {
			nextnodenr_ = 0;
			result = this;
		}
		return result;
	}

	void push_back(const TreeNode& node) {
		down_.emplace_back(new TreeNode(node));
		down_.back()->up = this;
	}

	size_t size() const { return down_.size(); }

	size_t Address() const { return address_; }

	void SetAddress(size_t addr) { address_ = addr; }

	const TensorDim& TDim() const { return tdim_; }

	void SetTDim(const TensorDim& tdim) { tdim_ = tdim; }

	const TreeNode& Down(size_t k) const {
		assert(k < down_.size());
		return *down_[k];
	}

	TreeNode& Down(size_t k) {
		assert(k < down_.size());
		return *down_[k];
	}

	void print(std::ostream& os = std::cout) const {
		content_.print(Layer(), os);
		for (const TreeNode *node : down_) {
			node->print(os);
		}
	}

	void GenInput(ostream& os = cout) const {
		// Indentation
		for (size_t i = 0; i < Layer() - 1; ++i) { os << "\t"; }
		if (IsLeaf()) {
			os << content_.n << "\t" << "6" << "\t" << content_.f << "\n";
		} else {
			os << content_.n << "\t-" << size() << "\n";
		}
		for (TreeNode *child : down_) {
			child->GenInput(os);
		}
	}

	bool IsLeaf() const { return (down_.empty()); }

	bool IsBottom() const { return (down_.size() == 1); }

	bool IsRoot() const { return (up == nullptr); }

	Path GetPath() const { return path_; }

	void SetPath(const Path& newpath) {
		path_ = newpath;
		for (size_t k = 0; k < down_.size(); ++k) {
			Path cpath = path_;
			cpath.emplace_back(k);
			TreeNode& child = Down(k);
			child.SetPath(cpath);
		}
	}

	void SetPhysicalNodes(size_t& next) {
		/// @TODO: Add a vector<size_t> with list of physical modes
		if (IsLeaf()) {
			content_.f = next;
			next++;
		} else {
			for (TreeNode *child : down_) {
				child->SetPhysicalNodes(next);
			}
		}
	}

	void SetPhysicalNodesScatter(size_t& nextl, size_t& nextr, bool& left) {
		if (IsLeaf()) {
			if (left) {
				content_.f = nextl;
				nextl++;
			} else {
				content_.f = nextr;
				nextr++;
			}
			left = !left;
		} else {
			for (TreeNode *child : down_) {
				child->SetPhysicalNodesScatter(nextl, nextr, left);
			}
		}
	}

	void UpdateTensorDim() {
		if( !IsLeaf()) {
			vector<size_t> ns;
			for (auto& child : down_) {
				ns.push_back(child->content_.n);
			}
			tdim_ = TensorDim(ns, content_.n);
			for (auto& child : down_) {
				child->UpdateTensorDim();
			}
		}
	}

	void MakeRoot() {
		up = nullptr;
		SetPath({0});
		UpdateTensorDim();
	}

	size_t Layer() const { return path_.size(); }

	NodeContent content_;
	TensorDim tdim_;
private:
	TreeNode *up;
	std::vector<TreeNode *> down_;

	int nextnodenr_;

	Path path_;

	size_t address_;
};


#endif //MCTDH_TREENODE_H
