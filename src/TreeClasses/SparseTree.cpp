//
// Created by Roman Ellerbrock on 2/2/20.
//
#include "TreeClasses/SparseTree.h"

SparseTree::SparseTree(const MLOcd& M, const Tree& tree, bool tail, bool inverse_tree)
	:SparseTree(M.targetLeaves(), tree, tail, inverse_tree) {}

SparseTree::SparseTree(const vector<size_t>& modes,
	const Tree& tree, bool tail, bool inverse_tree) {

	if (inverse_tree) {
		SparseTree stree(modes, tree, tail, false);
		initializeInverse(stree, tree);
	} else {
		SparseInitialize(modes, tree, tail);
	}
}

SparseTree::SparseTree(const SOPcd& sop, const Tree& tree, bool tail) {
	vector<size_t> actives;
	for (const MLOcd& M : sop) {
		for (size_t i = 0; i < M.size(); ++i) {
			size_t k = M.mode_(i);
			if (count(actives.begin(), actives.end(), k) == 0) {
				actives.push_back(k);
			}
		}
	}
	SparseInitialize(actives, tree, tail);
}

void SparseTree::SparseInitialize(const vector<size_t>& modes,
	const Tree& tree, bool tail) {

	/// Fill co_address with addresses in original TTBasis for every occuring node
	co_address_.clear();
	for (size_t k : modes) {
		const Leaf& phy = tree.GetLeaf(k);
		auto node = (const Node *) &phy.Up();
		const auto& beg = co_address_.begin();
		while (true) {
			size_t addr = node->Address();
			size_t count = co_address_.count(addr);
			// If node already there, increase its rank, otherwise add the node
			if (count == 0) {
				co_address_.insert(beg, pair<int, int>(node->Address(), 0));
			} else {
				co_address_[addr] += 1;
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
	nodes_.clear();
	for (auto& entry : co_address_) {
		entry.second = n++;
		const Node *node = &tree.GetNode(entry.first);
		nodes_.push_back(node);
	}

	if (!tail) {
		/// Go top-down and look for first node with more than one active children
		int m = nodes_.size() - 1;
		for (m = nodes_.size() - 1; m > 0; --m) {
			const Node& nodep = node(m);
			size_t NumActiveChildren = 0;
			for (size_t k = 0; k < nodep.nChildren(); ++k) {
				if (isActive(nodep.child(k))) { NumActiveChildren++; }
			}
			if (NumActiveChildren > 1) { break; }
		}
		if (m < (nodes_.size())) { m++; }
		/// Cut of tail
		for (size_t l = m; l < nodes_.size(); ++l) {
			const Node* node = nodes_[l];
			size_t addr = node->Address();
			co_address_.erase(addr);

		}
		nodes_ = vector<const Node *>(nodes_.begin(), nodes_.begin() + m);
	}
}

void SparseTree::print(const Tree& tree, ostream& os) const {
	for (const Node *node : *this) {
		node->info();
	}
}

void SparseTree::initializeInverse(const SparseTree& stree, const Tree& tree) {
	nodes_.clear();
	co_address_.clear();
	size_t addr = 0;
	for (const Node& node : tree) {
		if (!stree.isActive(node)) {
			nodes_.push_back(&node);
			co_address_[node.Address()] = addr;
			addr++;
		}
	}
}

