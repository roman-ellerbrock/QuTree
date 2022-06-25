//
// Created by Roman Ellerbrock on 6/22/22.
//

#ifndef SYMMATRIXTREE_H
#define SYMMATRIXTREE_H

using SymMatrixTree = pair<SparseMatrixTreecd, SparseMatrixTreecd>;
using SymMatrixTrees = vector<SymMatrixTree>;

/*class SymMatrixTree : public pair<SparseMatrixTreecd, SparseMatrixTreecd> {
public:
	explicit SymMatrixTree(const Tree& tree)
		: pair<SparseMatrixTreecd, SparseMatrixTreecd>(tree, tree) {
	}



};

using SymMatrixTrees = vector<SymMatrixTree>;
 */

#endif //SYMMATRIXTREE_H
