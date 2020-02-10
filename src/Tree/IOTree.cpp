//
// Created by Roman Ellerbrock on 2/10/20.
//

#include "IOTree.h"
#include "SpectralDecompositionTree.h"

namespace IOTree {

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const TTBasis& basis, ostream& os) {
		HoleMatrixTree<T> Rho(Psi, basis);
		SpectralDecompositionTree<T> specs(Psi, Rho, basis);
		specs.print(basis);
	}

	template <typename T>
	void Leafs(const TensorTree<T>& Psi, const TTBasis& basis, ostream& os) {

	}
}

template IOTree::Occupance<complex<double>>;
template IOTree::Occupance<double>;
