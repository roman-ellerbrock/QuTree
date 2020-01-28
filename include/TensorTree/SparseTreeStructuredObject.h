//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef MCTDH_SPARSETREESTRUCTUREDOBJECT_H
#define MCTDH_SPARSETREESTRUCTUREDOBJECT_H
#include "TreeMarker.h"
#include "TreeStructuredObject.h"

/**
 * \class SparseTreeStructuredObject
 *
 * \brief This is a central base class for objects that are sparsely tree-structured
 *
 * This is a sparse version of the TreeStructuredObject. It provides a completely
 * sparse handling for objects that are only present at some nodes.
 * */

template<class A>
class SparseTreeStructuredObject {
public:
	SparseTreeStructuredObject(const vector<size_t>& modes,
	const TTBasis& basis) : active(make_shared<TreeMarker>(modes, basis)) {
		Initialize(basis);
	}

	SparseTreeStructuredObject(shared_ptr<TreeMarker>& active_,
		const TTBasis& basis)
		: active(active_) {
		Initialize(basis);
	}

	virtual void Initialize(const TTBasis& basis) {
		attributes.resize(Active().size());
	}

	// Getter for the attribute to object a
	A& operator[](const Node& node) {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes.size());
		assert(Active(node));
		return attributes[address];
	}

	const A& operator[](const Node& node) const {
		size_t address = Active().SparseAddress(node);
		assert(address < attributes.size());
		assert(Active(node));
		return attributes[address];
	}

	typename vector<A>::const_iterator begin() const {
		return attributes.begin();
	}

	typename vector<A>::const_iterator end() const {
		return attributes.end();
	}

	size_t Size() const { return attributes.size(); }

	const TreeMarker& Active() const { return *active.get(); }

	size_t Active(const Node& node) const {
		return active->Active(node);
	}
	/**
	 * \brief Erase the data assigned to x.
	 * \param x Object which's attributes get deleted
	 *
	 * After erasing the attributes, the old structure needs to
	 * be updated (addresses become invalid).
	 */

protected:
	shared_ptr<TreeMarker> active;
	vector<A> attributes;
};

/*
SparseTreeStructuredObject(const MPO<T>& M,
	const TTBasis& basis)
	: active(make_shared<TreeMarker>(M, basis)) {
	Initialize(basis);
}
*/

#endif //MCTDH_SPARSETREESTRUCTUREDOBJECT_H
