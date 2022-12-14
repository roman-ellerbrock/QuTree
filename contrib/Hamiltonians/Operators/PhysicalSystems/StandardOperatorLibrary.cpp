//
// Created by Roman Ellerbrock on 6/25/20.
//

#include "StandardOperatorLibrary.h"

namespace Operator {

	SOPcd KineticEnergy(const Tree& tree) {
		/// Kinetic energy as Laplace operator
		LeafFuncd kin = &LeafInterface::applyKin;
		SOPcd H;
		for (size_t k = 0; k < tree.nLeaves(); ++k) {
			MLOcd M;
			M.push_back(kin, k);
			H.push_back(M, 1.);
		}
		return H;
	}

}

