//
// Created by Roman Ellerbrock on 2/10/20.
//

#include "IOTree.h"
#include "SpectralDecompositionTree.h"

namespace IOTree {

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os) {
		MatrixTree<T> Rho(tree);
		MatrixTreeFunctions::Contraction(Rho, Psi, tree, true);
		SpectralDecompositionTree<T> specs(Rho, tree);
		specs.print(tree);
	}

	template <typename T>
	FactorMatrix<T> LeafDensity(const TensorTree<T>& Psi, const MatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree) {
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
	void Leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os) {
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			auto rho_leaf = LeafDensity(Psi, Rho, leaf, tree);
			cout << "Leaf: " << l << "\n";
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

template void IOTree::Occupancy<complex<double>>(const TensorTree<complex<double>>& Psi, const Tree& tree, ostream& os);
template void IOTree::Occupancy<double>(const TensorTree<double>& Psi, const Tree& tree, ostream& os);

template void IOTree::Leafs<cd>(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Tree& tree, ostream& os);
template void IOTree::Leafs<d>(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Tree& tree, ostream& os);

template FactorMatrix<cd> IOTree::LeafDensity(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Leaf& leaf, const Tree& tree);
template FactorMatrix<d> IOTree::LeafDensity(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Leaf& leaf, const Tree& tree);


