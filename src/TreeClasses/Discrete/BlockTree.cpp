//
// Created by Roman Ellerbrock on 1/11/23.
//
#include "TreeClasses/Discrete/BlockTree.h"

size_t min(size_t a, size_t b) {
	return (a < b) ? a : b;
}

bool hasLabel(const Labels& L, const Label& lab) {
	return (std::find(L.begin(), L.end(), lab) != L.end());
}

ostream& operator<<(ostream& os, const LabelDimension& dim) {
	for (auto p : dim) {
		os << p.first << "\t" << p.second << endl;
	}
	return os;
}

Labels combine(const Labels& L, const Labels& R, const Range& range) {
	if (L.empty()) { return R; }
	if (R.empty()) { return L; }
	Labels res;
	for (size_t i = 0; i < L.size(); ++i) {
		for (size_t j = 0; j < R.size(); ++j) {
			Label new_label = L[i] + R[j];
			if (range.isAllowed(new_label)) {
				if (!hasLabel(res, new_label)) { res.push_back(new_label); }
			}
		}
	}
	std::sort(res.begin(), res.end());
	return res;
}

Labels label_combinations(const NodeAttribute<Labels>& up, const NodeAttribute<Labels>& down,
	const Range& range, const Node& node, const Node& hole) {

	Labels labels;
	if (!node.isBottomlayer()) {
		for (size_t i = 0; i < node.nChildren(); ++i) {
			const Node& child = node.child(i);
			if (child.address() == hole.address()) { continue; }
			labels = combine(labels, up[child], range);
		}
	}

	if (!node.isToplayer()) {
		if (node.address() != hole.address()) {
			labels = combine(labels, down[node], range);
		}
	}
	return labels;
}

size_t factorial(size_t n) {
	size_t res{1};
	for (size_t k = 1; k <= n; ++k) {
		res *= k;
	}
	return res;
}

size_t binomial(size_t k, size_t n) {
	return factorial(n) / (factorial(n - k) * factorial(k));
}

size_t labelDimension(const Label& label, size_t dimension, size_t n_hole) {
	size_t dim = binomial(label, n_hole);
	return min(dim, dimension);
}

