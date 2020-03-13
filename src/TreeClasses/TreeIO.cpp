//
// Created by Roman Ellerbrock on 2/10/20.
//

#include "TreeClasses/TreeIO.h"
#include "TreeClasses/SpectralDecompositionTree.h"

namespace TreeIO {

	template<typename T>
	void Output(const TensorTree<T>& Psi, const Tree& tree, ostream& os) {
		MatrixTree<T> Rho(tree);
		TreeFunctions::Contraction(Rho, Psi, tree, true);
		Occupancy(Psi, tree, os);
		Leafs(Psi, Rho, tree, os);
	}

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os) {
		MatrixTree<T> Rho(tree);
		TreeFunctions::Contraction(Rho, Psi, tree, true);
		SpectralDecompositionTree<T> specs(Rho, tree);
		specs.print(tree);
	}

	template <typename T>
	Matrix<T> LeafDensity(const TensorTree<T>& Psi, const MatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree) {
		const auto& node = (const Node&) leaf.Up();
		const auto& Phi = Psi[node];
		if (!node.isToplayer()) {
			const auto& rho = Rho[node];
			auto rhoPhi = multStateAB<T>(rho, Phi);
			return mHoleProduct(Phi, rhoPhi, 0);
		} else {
			return mHoleProduct(Phi, Phi, 0);
		}
	}

	template <typename T>
	void Leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os) {
		os << fixed;
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			auto rho_leaf = LeafDensity(Psi, Rho, leaf, tree);
			cout << "Leaf: " << l << "\n";
			double norm = abs(rho_leaf.Trace());
			for  (size_t i = 0; i < rho_leaf.Dim1(); ++i) {
				os << abs(rho_leaf(i, i)) / norm << "\t";
			}
			os << "\n";
		}
		os.flush();
		os << defaultfloat;
	}

}

typedef complex<double>  cd;
typedef double  d;

template void TreeIO::Occupancy<complex<double>>(const TensorTree<complex<double>>& Psi, const Tree& tree, ostream& os);
template void TreeIO::Occupancy<double>(const TensorTree<double>& Psi, const Tree& tree, ostream& os);

template void TreeIO::Output<d>(const TensorTree<d>& Psi, const Tree& tree, ostream& os);
template void TreeIO::Output<cd>(const TensorTree<cd>& Psi, const Tree& tree, ostream& os);

template void TreeIO::Leafs<cd>(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Tree& tree, ostream& os);
template void TreeIO::Leafs<d>(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Tree& tree, ostream& os);

template Matrix<cd> TreeIO::LeafDensity(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Leaf& leaf, const Tree& tree);
template Matrix<d> TreeIO::LeafDensity(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Leaf& leaf, const Tree& tree);


