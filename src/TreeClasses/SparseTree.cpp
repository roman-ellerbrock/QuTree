//
// Created by Roman Ellerbrock on 2/2/20.
//
#include "SparseTree.h"

void SparseTree::SparseInitialize(const vector<size_t>& modes,
	const Tree& tree, bool tail) {

	/// Fill co_address with addresses in original TTBasis for every occuring node
	co_address.clear();
	for (size_t k : modes) {
		const Leaf& phy = tree.GetLeaf(k);
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
			if (node->isToplayer()) {
				break;
			} else {
				node = &(node->parent());
			}
		}
	}

	/// Fill nodes vector with pointers to nodes and set sparse addresses
	size_t n = 0;
	nodes.clear();
	for (auto& entry : co_address) {
		entry.second = n++;
		const Node *node = &tree.GetNode(entry.first);
		nodes.push_back(node);
	}

	if (!tail) {
		/// Go top-down and look for first node with more than one active children
		for (int n = nodes.size() - 1; n > 0; --n) {
			const Node& node = MCTDHNode(n);
			size_t NumActiveChildren = 0;
			for (size_t k = 0; k < node.nChildren(); ++k) {
				if (Active(node.child(k))) { NumActiveChildren++; }
			}
			if (NumActiveChildren > 1) { break; }
		}
		n--;
		/// Cut of tail
		nodes = vector<const Node *>(nodes.begin(), nodes.begin() + n);
	}
}

void SparseTree::print(const Tree& tree, ostream& os) const {
	for (const Node *node : *this) {
		node->info();
	}
}
