//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef MCTDH_SPARSETREESTRUCTUREDOBJECT_H
#define MCTDH_SPARSETREESTRUCTUREDOBJECT_H
#include "TreeMarker.h"
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
	/// Construct object by providing a list of leaf-modes that are active.
	/// Will find path connecting the nodes.
	SparseTreeStructuredObject(const vector<size_t>& modes,
		const TTBasis& basis)
		: active(make_shared<TreeMarker>(modes, basis)) {
		SparseTreeStructuredObject::Initialize(basis);
	}

	/// Construct obejct for previously marked active nodes
	SparseTreeStructuredObject(shared_ptr<TreeMarker>& active_,
		const TTBasis& basis)
		: active(active_) {
		SparseTreeStructuredObject::Initialize(basis);
	}

	/// Allocate memory. Call only after initializing TreeMarker. Requires default constructor.
	virtual void Initialize(const TTBasis& basis) {
		attributes.resize(Active().size());
	}

	/// Getter for the attribute at an active node. Crashes if called for inactive nodes.
	A& operator[](const Node& node) {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes.size());
		assert(Active(node));
		return attributes[address];
	}

	/// Getter for the attribute at an active node. Crashes if called for inactive nodes.
	const A& operator[](const Node& node) const {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes.size());
		assert(Active(node));
		return attributes[address];
	}

	/// Swipe bottom-up over attributes at every active node
	typename vector<A>::const_iterator begin() const {
		return attributes.begin();
	}

	/// Swipe bottom-up over attributes at every active node
	typename vector<A>::const_iterator end() const {
		return attributes.end();
	}

	/// Number of active nodes
	size_t Size() const { return attributes.size(); }

	/// Getter to TreeMarker
	const TreeMarker& Active() const { return *active.get(); }

	/// Check whether node is active
	size_t Active(const Node& node) const {
		return active->Active(node);
	}

protected:
	/// TreeMarker which marks whether a node is active or not
	shared_ptr<TreeMarker> active;
	/// Attributes only at every active node
	vector<A> attributes;
};

#endif //MCTDH_SPARSETREESTRUCTUREDOBJECT_H
