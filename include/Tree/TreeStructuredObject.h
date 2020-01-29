#pragma once
#include "TensorTreeBasis/Node.h"

/**
 * \class TreeStructuredObject
 * \brief This is a key structure-giving class.
 * Each element of a set of elements A is connected to one
 * element of a set of elements B.
 *
 * Assume that there is a set X with objects {x_1,...,x_N}.
 * A AttributiveSet is a set {a_1,...,a_N} where each
 * object a_i is asigned to an object x_i.
 * */

template<class A>
class TreeStructuredObject {
public:
	// Getter for the attribute to object a
	A& operator[](const Node& x) {
		size_t address = x.Address();
		assert(address < attributes.size());
		return attributes[address];
	}

	const A& operator[](const Node& x) const {
		size_t address = x.Address();
		assert(address < attributes.size());
		return attributes[address];
	}

	typename vector<A>::iterator begin() {
		return attributes.begin();
	}

	typename vector<A>::iterator end() {
		return attributes.end();
	}

	typename vector<A>::const_iterator begin() const {
		return attributes.begin();
	}

	typename vector<A>::const_iterator end() const {
		return attributes.end();
	}

	size_t size() const { return attributes.size(); }

	/**
	 * \brief Erase the data assigned to x.
	 * \param x Object which's attributes get deleted
	 *
	 * After erasing the attributes, the old structure needs to
	 * be updated (addresses become invalid).
	 */
	void erase(const Node& x) {
		attributes.erase(attributes.begin() + x.Address());
	}

protected:
	// Ordered list (vector) of Attributes a_i
	vector<A> attributes;
};


