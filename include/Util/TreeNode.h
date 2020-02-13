//
// Created by Roman Ellerbrock on 2019-04-22.
//

#ifndef MCTDH_TREENODE_H
#define MCTDH_TREENODE_H
#include <utility>
#include <vector>
#include <iostream>
#include "Core/TensorDim.h"

typedef vector<size_t> Path;

void print(const Path& p, ostream& os = cout) {
		os << "{ ";
	if (p.empty()) {
		os << "{ }" << endl;
		return;
	}
	os << p[0];
	for (size_t i = 1; i < p.size(); ++i) {
		os << ", " << p[i];
	}
	os << " }" << endl;
}

class NodeContent {
public:
	NodeContent()
		: f(0), n(0), nextnodenr_(0), address_(0) {}

	explicit NodeContent(size_t n_, size_t f_ = 1)
		: f(f_), n(n_), nextnodenr_(0), address_(0) {}

	~NodeContent() = default;

	void print(size_t indent = 0, std::ostream& os = std::cout) const {
		for (size_t i = 0; i < indent; ++i) {
			os << "\t";
		}
		os << f << " " << n << std::endl;
	}

	size_t f, n;

	int nextnodenr_;

	size_t address_;

	TensorDim tdim_;

	Path path_;
};

class TreeNode {
public:

	/**
	 * A simple tree-node class.
	 */

	/// Rule of five-section (constructors & destructors)
	TreeNode()
		: up_(nullptr) {}

	explicit TreeNode(size_t n)
		: up_(nullptr), content_(n) {}

	explicit TreeNode(NodeContent con)
		: up_(nullptr), content_(con) {}

	~TreeNode() = default;

	/// Copy constructor
	TreeNode(const TreeNode& node)
		: up_(node.up_), content_(node.content_) {

		for (const TreeNode *child : node.down_) {
			down_.push_back(new TreeNode(*child));
		}
		for (TreeNode *child : down_) {
			child->up_ = this;
		}
	}

	/// Move constructor
	TreeNode(TreeNode&& node) noexcept
		: up_(node.up_), down_(move(node.down_)), content_(node.content_) {

		for (TreeNode *child : down_) {
			child->up_ = this;
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

		up_ = old.up_;
		down_ = move(old.down_);
		content_ = old.content_;

		for (TreeNode *child : down_) {
			child->up_ = this;
		}
		return *this;
	}

	/// Member functions
	TreeNode *NextNode() {
		TreeNode *result;
		if (content_.nextnodenr_ < down_.size()) {
			result = down_[content_.nextnodenr_]->NextNode();
			if (result == down_[content_.nextnodenr_]) {
				++content_.nextnodenr_;
			}
		} else {
			content_.nextnodenr_ = 0;
			result = this;
		}
		return result;
	}

	void push_back(const TreeNode& node) {
		down_.emplace_back(new TreeNode(node));
		down_.back()->up_ = this;
	}

	size_t size() const { return down_.size(); }

	size_t Address() const { return content_.address_; }

	void SetAddress(size_t addr) { content_.address_ = addr; }

	const TensorDim& TDim() const { return content_.tdim_; }

	void SetTDim(const TensorDim& tdim) { content_.tdim_ = tdim; }

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

	bool IsLeaf() const { return (down_.empty()); }

	bool IsBottom() const { return (down_.size() == 1); }

	bool IsRoot() const { return (up_ == nullptr); }

	size_t Layer() const { return content_.path_.size(); }

	Path GetPath() const { return content_.path_; }

	void SetPath(const Path& newpath) {
		content_.path_ = newpath;
		for (size_t k = 0; k < down_.size(); ++k) {
			Path cpath = content_.path_;
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

	void UpdateTensorDim() {
		if( !IsLeaf()) {
			vector<size_t> ns;
			for (auto& child : down_) {
				ns.push_back(child->content_.n);
			}
			content_.tdim_ = TensorDim(ns, content_.n);
			for (auto& child : down_) {
				child->UpdateTensorDim();
			}
		}
	}

	void MakeRoot() {
		up_ = nullptr;
		SetPath({0});
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

	NodeContent content_;
private:
	TreeNode *up_;
	std::vector<TreeNode *> down_;

};


#endif //MCTDH_TREENODE_H
