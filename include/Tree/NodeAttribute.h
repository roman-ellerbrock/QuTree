#pragma once
#include "TreeHandling/Node.h"

/**
 * \defgroup Tree
 * \brief Group for classes related to tree handling.
 *
 * Classes in this group use a TTBasis to generate new classes
 * that use and work with the tree structure provided by TTBasis.
 * The classes in this group should either inherit from
 * TreeStructureObject (for classes that have some object at
 * EVERY node in a tree) or SparseTreeStructuredObject (for classes
 * that have some object at a SUBSET of all nodes in the tree).
 */

template<class A>
class NodeAttribute
/**
 * \class NodeAttribute
 * \ingroup Tree
 * \brief Base class for creating objects at every node in a tree.
 *
 * For every node in a TTBasis, there is one corresponding Object
 * of class A. Inherit from this class for fast prototyping.
 * Make sure that attributes_ gets cleared and filled in a constructor.
 * Automatically provides iterators and bracket operators for Node objects.
 * */
{
public:
	// Getter for the attribute to object a
	A& operator[](const Node& x) {
		size_t address = x.Address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	const A& operator[](const Node& x) const {
		size_t address = x.Address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	typename vector<A>::iterator begin() {
		return attributes_.begin();
	}

	typename vector<A>::iterator end() {
		return attributes_.end();
	}

	typename vector<A>::const_iterator begin() const {
		return attributes_.begin();
	}

	typename vector<A>::const_iterator end() const {
		return attributes_.end();
	}

	size_t size() const { return attributes_.size(); }

	/**
	 * \brief Erase the data assigned to x_.
	 * \param x Object which's attributes_ get deleted
	 *
	 * After erasing the attributes_, the old structure needs to
	 * be updated (addresses become invalid).
	 */
	void erase(const Node& x) {
		attributes_.erase(attributes_.begin() + x.Address());
	}

protected:
	// Ordered list (vector) of Attributes a_i
	vector<A> attributes_;
};


