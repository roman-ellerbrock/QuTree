//
// Created by Roman Ellerbrock on 2/10/20.
//

#include "IOTree.h"
#include "SpectralDecompositionTree.h"

namespace IOTree {

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const TTBasis& basis, ostream& os) {
		HoleMatrixTree<T> Rho(Psi, basis);
		SpectralDecompositionTree<T> specs(Rho, basis);
		specs.print(basis);
	}

	template <typename T>
	FactorMatrix<T> LeafDensity(const TensorTree<T>& Psi, const HoleMatrixTree<T>& Rho,
		const Leaf& leaf, const TTBasis& basis) {
		const auto& node = (const Node&) leaf.Up();
		const auto& Phi = Psi[node];
		if (!node.IsToplayer()) {
			const auto& rho = Rho[node];
			auto rhoPhi = multStateAB<T>(rho, Phi);
			return HoleProduct(Phi, rhoPhi, 0);
		} else {
			return HoleProduct(Phi, Phi, 0);
		}
	}

	template <typename T>
	void Leafs(const TensorTree<T>& Psi, const HoleMatrixTree<T>& Rho, const TTBasis& basis, ostream& os) {
		for (size_t l = 0; l < basis.nLeaves(); ++l) {
			const Leaf& leaf = basis.GetLeaf(l);
			auto rho_leaf = LeafDensity(Psi, Rho, leaf, basis);
			cout << "Leaf " << l << endl;
			for  (size_t i = 0; i < rho_leaf.Dim(); ++i) {
				os << abs(rho_leaf(i, i)) << "\t";
			}
			os << "\n";
		}
		os.flush();
	}
}

typedef complex<double>  cd;
typedef double  d;

template void IOTree::Occupancy<complex<double>>(const TensorTree<complex<double>>& Psi, const TTBasis& basis, ostream& os);
template void IOTree::Occupancy<double>(const TensorTree<double>& Psi, const TTBasis& basis, ostream& os);

template void IOTree::Leafs<cd>(const TensorTree<cd>& Psi, const HoleMatrixTree<cd>& Rho, const TTBasis& basis, ostream& os);
template void IOTree::Leafs<d>(const TensorTree<d>& Psi, const HoleMatrixTree<d>& Rho, const TTBasis& basis, ostream& os);

template FactorMatrix<cd> IOTree::LeafDensity(const TensorTree<cd>& Psi, const HoleMatrixTree<cd>& Rho, const Leaf& leaf, const TTBasis& basis);
template FactorMatrix<d> IOTree::LeafDensity(const TensorTree<d>& Psi, const HoleMatrixTree<d>& Rho, const Leaf& leaf, const TTBasis& basis);


