//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef SPARSETREESTRUCTUREDOBJECT_H
#define SPARSETREESTRUCTUREDOBJECT_H
#include "SubTree.h"
#include "TreeStructuredObject.h"


template<class A>
class SparseTreeStructuredObject
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
	SparseTreeStructuredObject(const vector<size_t>& modes,
		const TTBasis& tree)
		: active_(make_shared<SubTree>(modes, tree)) {
		SparseTreeStructuredObject::Initialize(tree);
	}

	/// Construct obejct for previously marked active_ nodes
	SparseTreeStructuredObject(shared_ptr<SubTree>& active_,
		const TTBasis& tree)
		: active_(active_) {
		SparseTreeStructuredObject::Initialize(tree);
	}

	/// Allocate memory. Call only after initializing TreeMarker. Requires default constructor.
	virtual void Initialize(const TTBasis& tree) {
		attributes_.resize(Active().size());
	}

	/// Getter for the attribute at an active_ node. Crashes if called for inactive nodes.
	A& operator[](const Node& node) {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes_.size());
		assert(Active(node));
		return attributes_[address];
	}

	/// Getter for the attribute at an active_ node. Crashes if called for inactive nodes.
	const A& operator[](const Node& node) const {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes_.size());
		assert(Active(node));
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
	size_t Size() const { return attributes_.size(); }

	/// Getter to TreeMarker
	const SubTree& Active() const { return *active_.get(); }

	/// Check whether node is active_
	size_t Active(const Node& node) const {
		return active_->Active(node);
	}

protected:
	/// TreeMarker which marks whether a node is active_ or not
	shared_ptr<SubTree> active_;
	/// Attributes only at every active_ node
	vector<A> attributes_;
};

#endif //SPARSETREESTRUCTUREDOBJECT_H
