//
//

#include "DVR/SymXMatrixTrees.h"

SOPcd symXsop(const Tree& tree) {
	LeafFuncd x = &LeafInterface::applyX;
	LeafFuncd I = &LeafInterface::identity;
	SOPcd xops;
	for (size_t l = 0; l < tree.nLeaves(); ++l) {
		const Leaf& leaf = tree.getLeaf(l);
		size_t mode = leaf.mode();
		MLOcd M(x, mode);
//		for (size_t i = 0; i < tree.nLeaves(); ++i) {
//			M.push_back(I, i);
//		}
		xops.push_back(M, 1.);
	}
	return xops;
}

