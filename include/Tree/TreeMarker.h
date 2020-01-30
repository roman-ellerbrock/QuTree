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

/**
 * \class TreeMarker
 * \brief This class marks active upper nodes, given a list of physical modes
 *
 * */

class TreeMarker {
public:

	TreeMarker(const vector<size_t>& modes,
		const TTBasis& basis) {
		SparseInitialize(modes, basis);
	}

	void SparseInitialize(const vector<size_t>& modes, const TTBasis& basis) {
		co_address.clear();
		for (size_t k : modes) {
			const Leaf& phy = basis.GetLeaf(k);
			auto node = (const Node *) &phy.Up();
			const auto& beg = co_address.begin();
			while (true) {
				size_t addr = node->Address();
				size_t count = co_address.count(addr);
				// If node already there, increase its rank, otherwise add the node
				if (count == 0) {
					co_address.insert(beg, pair<int, int>(node->Address(), 0));
				} else {
					co_address[addr] += 1;
				}
				if (node->IsToplayer()) {
					break;
				} else {
					node = &(node->Up());
				}
			}
		}

		size_t n = 0;
		nodes.clear();
		for (auto& entry : co_address) {
			entry.second = n++;
			const Node *node = &basis.GetNode(entry.first);
			nodes.push_back(node);
		}
	}

	size_t Active(const Node& node) const {
		size_t count = co_address.count(node.Address());
//		return !(count == 0); ?
		if (0 == count) {
			return false;
		} else {
			return true;
		}
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

	void Initialize(const TTBasis& basis) {
//        attributes.resize(basis.nTotalNodes());
	}

	void print(const TTBasis& basis, ostream& os = cout) const{
		for (const Node* node : *this) {
			node->info();
		}
	}

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
