//
// Created by Roman Ellerbrock on 10/21/21.
//

#ifndef OVERLAPS_H
#define OVERLAPS_H
#include "TreeClasses/TensorTree.h"

enum OverlapType {
	full,
	diagonal
};

void wavefunctionOverlap(const string& file1, const string& file2, const Tree& tree, OverlapType = diagonal);

void wavefunctionOverlap(const vector<TensorTreecd>& Psi,
	const vector<TensorTreecd>& Chi, const Tree& tree, OverlapType = full);

#endif //OVERLAPS_H
