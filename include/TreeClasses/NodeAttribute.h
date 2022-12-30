#pragma once
#include "TreeShape/Node.h"
#include "TreeShape/Edge.h"
#include "TreeShape/Tree.h"

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
		size_t address = x.address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	const A& operator[](const Node& x) const {
		size_t address = x.address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	A& operator[](const Edge& e) {
		const Node& x = e.down();
		size_t address = x.address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	const A& operator[](const Edge& e) const {
		const Node& x = e.down();
		size_t address = x.address();
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

	void clear() { attributes_.clear(); }

	void emplace_back(A&& a) {
		attributes_.template emplace_back(a);
	}

	/**
	 * \brief Erase the data assigned to x_.
	 * \param x Object which's attributes_ get deleted
	 *
	 * after erasing the attributes_, the old structure needs to
	 * be updated (addresses become invalid).
	 */
	void erase(const Node& x) {
		attributes_.erase(attributes_.begin() + x.address());
	}

protected:
	// Ordered list (vector) of Attributes a_i
	vector<A> attributes_;
};

template <class A>
void print(const NodeAttribute<A>& ATree, const Tree& tree) {
	for (size_t i = 0; i < ATree.size(); ++i) {
		const auto& node_ = tree.begin() + i;
		const Node& node = *node_;
		node.info();
		ATree[node].print();
	}
}
