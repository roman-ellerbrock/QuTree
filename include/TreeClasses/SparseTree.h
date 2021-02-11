//
// Created by Roman Ellerbrock on 2019-04-26.
//

#ifndef MCTDH_TREEMARKER_H
#define MCTDH_TREEMARKER_H
#include "NodeAttribute.h"
#include "TreeShape/Tree.h"
#include "TreeOperators/MultiLeafOperator.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include <map>
#include <chrono>


class SparseTree
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
	SparseTree() = default;

	SparseTree(const vector<size_t>& modes,
		const Tree& tree, bool tail = true, bool inverse_tree = false);

	SparseTree(const MLOcd& M, const Tree& tree, bool tail = true, bool inverse_tree = false);

	SparseTree(const SOPcd& sop, const Tree& tree, bool tail = true);

	void SparseInitialize(const vector<size_t>& modes,
		const Tree& tree, bool tail = true);

	void initializeInverse(const SparseTree& stree,
		const Tree& tree);

	size_t isActive(const Node& node) const {
		size_t count = co_address_.count(node.Address());
		return (count != 0);
	}

	size_t size() const { return nodes_.size(); }

	vector<const Node *>::const_iterator begin() const {
		return nodes_.begin();
	}

	vector<const Node *>::const_iterator end() const {
		return nodes_.end();
	}

	const Node& node(size_t n) const {
		assert(n < nodes_.size());
		return *nodes_[n];
	}

	size_t sparseAddress(const Node& node) const {
		size_t addr = node.Address();
		return co_address_.at(addr);
	}

	void print(const Tree& tree, ostream& os = cout) const;

protected:
	vector<const Node *> nodes_;
	map<size_t, size_t> co_address_;
};

/*
TreeMarker(const MultiLeafOperator<T>& M,
	const TTBasis& tree) {
	vector<size_t> modes;
	for (size_t k = 0; k < M.size(); ++k) {
		modes.push_back(M.Mode(k));
	}
	SparseInitialize(modes, tree);
}

TreeMarker(const SOP& sop, const TTBasis& tree) {
	vector<size_t> actives;
	for (const MultiLeafOperator<T>& M : sop) {
		for (size_t i = 0; i < M.size(); ++i) {
			size_t k = M.Mode(i);
			if (!count(actives.begin(), actives.end(), k)) {
				actives.push_back(k);
			}
		}
	}
	SparseInitialize(actives, tree);
}
*/

#endif //MCTDH_TREEMARKER_H
