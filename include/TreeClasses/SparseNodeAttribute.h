//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef SPARSENODEATTRIBUTE_H
#define SPARSENODEATTRIBUTE_H
#include "SparseTree.h"
#include "NodeAttribute.h"

template<class A>
class SparseNodeAttribute
/**
 * \class SparseTreeStructuredObject
 * \ingroup Tree
 * \brief This is a central base class for objects that are sparsely tree-structured
 *
 * This is a sparse version of the TreeStructuredObject. It provides a
 * sparse handling for objects that are only present at some nodes.
 * */
{
public:
	/// Construct object by providing a list of leaf-modes that are active_.
	/// Will find path connecting the nodes.
	SparseNodeAttribute(const vector<size_t>& modes,
		const Tree& tree, bool tail = true, bool inverse = false)
		: active_(make_shared<SparseTree>(modes, tree, tail, inverse)) {
		SparseNodeAttribute::initialize(tree);
	}

	/// Construct object for previously marked active_ nodes
	SparseNodeAttribute(shared_ptr<SparseTree> active_,
		const Tree& tree)
		: active_(move(active_)) {
		SparseNodeAttribute::initialize(tree);
	}

	/// Construct object for previously marked active_ nodes
	SparseNodeAttribute(const SparseTree& stree, const Tree& tree)
		: active_(make_shared<SparseTree>(stree)) {
		SparseNodeAttribute::initialize(tree);
	}

	/// Construct object for previously marked active_ nodes
	explicit SparseNodeAttribute(const Tree& tree)
		: active_(make_shared<SparseTree>(tree)) {
		SparseNodeAttribute::initialize(tree);
	}

	/// Allocate memory. Call only after initializing TreeMarker. Requires default constructor.
	virtual void initialize(const Tree& tree) {
		attributes_.resize(sparseTree().size());
	}

	/// Getter for the attribute at an active_ node. Crashes if called for inactive nodes.
	A& operator[](const Node& node) {
		size_t address = sparseTree().sparseAddress(node);
		assert(address < attributes_.size());
		assert(isActive(node));
		return attributes_[address];
	}

	/// Getter for the attribute at an active_ node. Crashes if called for inactive nodes.
	const A& operator[](const Node& node) const {
		size_t address = sparseTree().sparseAddress(node);
		assert(address < attributes_.size());
		assert(isActive(node));
		return attributes_[address];
	}

	/// Swipe bottom-up over attributes_ at every active_ node
	typename vector<A>::const_iterator begin() const {
		return attributes_.begin();
	}

	/// Swipe bottom-up over attributes_ at every active_ node
	typename vector<A>::const_iterator end() const {
		return attributes_.end();
	}

	/// Number of active_ nodes
	[[nodiscard]] size_t size() const { return attributes_.size(); }

	/// Getter to TreeMarker
	[[nodiscard]] const SparseTree& sparseTree() const { return *active_; }

	/// Check whether node is active_
	[[nodiscard]] size_t isActive(const Node& node) const {
		return active_->isActive(node);
	}

protected:
	/// TreeMarker which marks whether a node is active_ or not
	shared_ptr<SparseTree> active_;
	/// Attributes only at every active_ node
	vector<A> attributes_;
};

#endif //SPARSENODEATTRIBUTE_H
