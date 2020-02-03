//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef MCTDH_TREEMARKER_H
#define MCTDH_TREEMARKER_H
#include "TreeStructuredObject.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "MultiParticleOperator.h"
#include "SumOfProductsOperator.h"
#include <map>
#include <chrono>


class TreeMarker
/**
 * \class TreeMarker
 * \brief This class marks a subset of active Nodes in a tree.
 *
 * The class is used to mark Nodes when working with sparseness in
 * tree structure. Typically, nodes are marked by providing a list
 * of active leaves. In this case, the TreeMarker searches for the
 * Nodes connectinb the leaves and saving the corresponding Node
 * pointers in a list.
 * co_address stores the mapping of the global Node address in
 * TTBasis to the sparse address.
 * */
{
public:

	TreeMarker(const vector<size_t>& modes,
		const TTBasis& basis, bool tail = true) {
		SparseInitialize(modes, basis, tail);
	}

	void SparseInitialize(const vector<size_t>& modes,
		const TTBasis& basis, bool tail = true);

	size_t Active(const Node& node) const {
		size_t count = co_address.count(node.Address());
		return (count != 0);
	}

	size_t size() const { return nodes.size(); }

	vector<const Node *>::const_iterator begin() const {
		return nodes.begin();
	}

	vector<const Node *>::const_iterator end() const {
		return nodes.end();
	}

	const Node& MCTDHNode(size_t n) const {
		assert(n < nodes.size());
		return *nodes[n];
	}

	size_t SparseAddress(const Node& node) const {
		size_t addr = node.Address();
		return co_address.at(addr);
	}

	void print(const TTBasis& basis, ostream& os = cout) const;

protected:
	vector<const Node *> nodes;
	map<size_t, size_t> co_address;
};

/*
TreeMarker(const MultiParticleOperator<T>& M,
	const TTBasis& basis) {
	vector<size_t> modes;
	for (size_t k = 0; k < M.size(); ++k) {
		modes.push_back(M.Mode(k));
	}
	SparseInitialize(modes, basis);
}

TreeMarker(const SOP& sop, const TTBasis& basis) {
	vector<size_t> actives;
	for (const MultiParticleOperator<T>& M : sop) {
		for (size_t i = 0; i < M.size(); ++i) {
			size_t k = M.Mode(i);
			if (!count(actives.begin(), actives.end(), k)) {
				actives.push_back(k);
			}
		}
	}
	SparseInitialize(actives, basis);
}
*/

#endif //MCTDH_TREEMARKER_H
