//
// Created by Roman Ellerbrock on 1/11/23.
//
#include "TreeClasses/Discrete/BlockTree.h"

bool hasLabel(const Labels& L, const Label& lab) {
	return (std::find(L.begin(), L.end(), lab) != L.end());
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

Labels label_combinations(const LabelTree& up, const LabelTree& down,
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

