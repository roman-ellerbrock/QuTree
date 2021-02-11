//
// Created by Roman Ellerbrock on 2/10/20.
//

#include "TreeClasses/TreeIO.h"
#include "TreeClasses/SpectralDecompositionTree.h"

namespace TreeIO {
	void Status(size_t it, size_t max, size_t freq, size_t length) {
		double ptg = it * 1. / ((double) max);
		auto num = (size_t) (ptg * length);
		if (!(it % freq)) {
			string msg = "Progress : [";
			for (size_t i = 0; i < num; ++i) {
				msg.append("#");
			}
			for (size_t i = num; i < length - 1; ++i) {
				msg.append(" ");
			}
			msg.append("]");
			cout << "\r" << msg;
			cout.flush();
		}
	}

	void StatusTime(size_t it, size_t max, size_t freq, size_t length,
		chrono::high_resolution_clock::time_point& t1, chrono::high_resolution_clock::time_point& t2,
		chrono::microseconds& time) {

		double ptg = it * 1. / ((double) max);
		auto num = (size_t) (ptg * length);
		if (!(it % freq)) {

			string clear
				("                                                                                                     ");
			cout << "\r" << clear;
			cout.flush();
			t2 = chrono::high_resolution_clock::now();
			time += chrono::duration_cast<chrono::microseconds>(t2 - t1);
			t1 = t2;
			double ttot = time.count() / 1000000.;
			string msg = "Total time: " + to_string(time.count() / 1000000.) + "s\t";
			msg.append("Time per operation: " + to_string(time.count() / (it * 1000.)) + "ms\t");
			msg.append("Progress : [");
			for (size_t i = 0; i < num; ++i) {
				msg.append("#");
			}
			for (size_t i = num; i < length - 1; ++i) {
				msg.append(" ");
			}
			msg.append("]");
			double project = ttot * length / (1. * num);
			cout << " | " << project;
			cout << "\r" << msg;
			cout.flush();
		}
	}

	void expectationValues(const TensorTreecd& Psi, const MatrixTreecd& rho,
		const Tree& tree, ostream& os) {

		// This routine calculates expectation values for the bottom-layer_ using
		// single-particle density matrices
		// <h>^k=tr(rho_k h)/tr(rho_k)

		// x,p-Operator
		function<void(const LeafInterface&, Tensorcd&, const Tensorcd&)>
			x = &LeafInterface::applyX;
		function<void(const LeafInterface&, Tensorcd&, const Tensorcd&)>
			p = &LeafInterface::applyP;

		// Calculate expectation values for each physica mode
		for (int k = 0; k < tree.nLeaves(); k++) {
			// References to the Physical coordinate and bottom-layer_ node
			const Leaf& phys = tree.GetLeaf(k);
			const LeafInterface& grid = phys.PrimitiveGrid();
			auto& node = (Node&) phys.Up();
			cout << "mode=" << phys.Mode() << endl;

			// Get SPF
			Tensorcd Acoeff = Psi[node];

			// state-multiply with single-particle density matrix
			const Matrixcd& dens = rho[node];
			double norm = abs(real(dens.trace()));
			Tensorcd hAcoeff = multStateAB(dens, Acoeff);

			// Calculate <x>
			{
				Tensorcd xhAcoeff(Acoeff.shape());
				Tensorcd xxhAcoeff(Acoeff.shape());
				x(grid, xhAcoeff, hAcoeff);
				x(grid, xxhAcoeff, xhAcoeff);

				// tr(rho*x)
				Matrixcd xmat = Acoeff.dotProduct(xhAcoeff);
				double xexp = real(xmat.trace() / norm);
				Matrixcd x2mat = Acoeff.dotProduct(xxhAcoeff);
				double x2exp = real(x2mat.trace() / norm);
				double dx = sqrt(abs(x2exp - xexp * xexp));
				// tr(rho)
				cout << "<x>=" << xexp << "    ";
				cout << "<dx>=" << dx << "    " << endl;
			}
			// occupation
			if (phys.Type() != 1) {
				Tensorcd nAcoeff(hAcoeff);
				Tensorcd n2Acoeff(hAcoeff);
				for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
					for (int i = 0; i < Acoeff.shape().lastBefore(); i++) {
						nAcoeff(i, n) *= i;
						n2Acoeff(i, n) *= i * i;
					}
				}
				Matrixcd nmat = Acoeff.dotProduct(nAcoeff);
				double nexp = real(nmat.trace() / norm);
				Matrixcd n2mat = Acoeff.dotProduct(n2Acoeff);
				double n2exp = real(n2mat.trace() / norm);
				double dn = sqrt(abs(n2exp - nexp * nexp));
//			cout << "<n>=" << nexp << "    <dn>=" << dn << endl;
			}
			// Calculate <p>
			if (phys.Type() != 2) {
				Tensorcd phAcoeff(hAcoeff.shape());
				Tensorcd pphAcoeff(hAcoeff.shape());
				p(grid, phAcoeff, hAcoeff);
				p(grid, pphAcoeff, phAcoeff);
				// tr(rho*x)
				Matrixcd pmat = Acoeff.dotProduct(phAcoeff);
				double pexp = real(pmat.trace() / norm);
				Matrixcd p2mat = Acoeff.dotProduct(pphAcoeff);
				double p2exp = real(p2mat.trace() / norm);
				double dp = sqrt(abs(p2exp - pexp * pexp));
				// tr(rho)
				cout << "<p>=" << pexp << "    ";
				cout << "<dp>=" << dp << endl;
			}
		}
	}

	void Output(const TensorTreecd& Psi, const Tree& tree, ostream& os) {
		MatrixTreecd Rho(tree);
		TreeFunctions::contraction(Rho, Psi, tree, true);
		Occupancy(Psi, tree, os);
		expectationValues(Psi, Rho, tree, os);
		Leafs(Psi, Rho, tree, os);
	}

	template<typename T>
	void Occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os) {
		MatrixTree<T> Rho(tree);
		TreeFunctions::contraction(Rho, Psi, tree, true);
		SpectralDecompositionTree<T> specs(Rho, tree);
		specs.print(tree);
	}

	template<typename T>
	Matrix<T> LeafDensity(const TensorTree<T>& Psi, const MatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree) {
		const auto& node = (const Node&) leaf.Up();
		const auto& Phi = Psi[node];
		if (!node.isToplayer()) {
			const auto& rho = Rho[node];
			auto rhoPhi = multStateAB<T>(rho, Phi);
			return contraction(Phi, rhoPhi, 0);
		} else {
			return contraction(Phi, Phi, 0);
		}
	}

	template<typename T>
	Matrix<T> LeafDensity(const TensorTree<T>& Psi, const SparseMatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree) {
		const auto& node = (const Node&) leaf.Up();
		const SparseTree& stree = Rho.sparseTree();
		assert(stree.Active(node));
		const auto& Phi = Psi[node];
		if (!node.isToplayer()) {
			const auto& rho = Rho[node];
			auto rhoPhi = multStateAB<T>(rho, Phi);
			return contraction(Phi, rhoPhi, 0);
		} else {
			return contraction(Phi, Phi, 0);
		}
	}

	template<typename T>
	void Leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os) {
		os << fixed;
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			auto rho_leaf = LeafDensity(Psi, Rho, leaf, tree);
			cout << "Leaf: " << l << "\n";
			double norm = abs(rho_leaf.trace());
			for (size_t i = 0; i < rho_leaf.dim1(); ++i) {
				os << abs(rho_leaf(i, i)) / norm << "\t";
				if ((i + 1) % 8 == 0) { os << "\n"; }
			}
			os << "\n";
		}
		os.flush();
		os << defaultfloat;
	}

	template <typename T>
	void EntropyMap(const TensorTree<T>& Psi, const Tree& tree) {
		auto rho = TreeFunctions::contraction(Psi, tree, true);
		SpectralDecompositionTree<T> X = diagonalize(rho);
		for (const SpectralDecomposition<T>& x : X) {
			const auto& occ = x.second;
			double s = 0.;
			for (size_t i = 0; i < occ.dim(); ++i) {
				s -= occ(i) * log(occ(i));
			}
		}
	}

	template<class A>
	void print(const vector<A>& vec) {
		for (const auto& element : vec) {
			element.print();
		}
	}
}

typedef complex<double> cd;

typedef double d;

template void TreeIO::Occupancy<complex<double>>(const TensorTree<complex<double>>& Psi, const Tree& tree, ostream& os);
template void TreeIO::Occupancy<double>(const TensorTree<double>& Psi, const Tree& tree, ostream& os);

//template void TreeIO::Output<d>(const TensorTree<d>& Psi, const Tree& tree, ostream& os);
//template void TreeIO::Output(const TensorTree<cd>& Psi, const Tree& tree, ostream& os);

template void TreeIO::Leafs<cd>(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Tree& tree, ostream& os);
template void TreeIO::Leafs<d>(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Tree& tree, ostream& os);

template Matrix<cd>
TreeIO::LeafDensity(const TensorTree<cd>& Psi, const MatrixTree<cd>& Rho, const Leaf& leaf, const Tree& tree);
template Matrix<d>
TreeIO::LeafDensity(const TensorTree<d>& Psi, const MatrixTree<d>& Rho, const Leaf& leaf, const Tree& tree);

template Matrix<cd>
TreeIO::LeafDensity(const TensorTree<cd>& Psi, const SparseMatrixTree<cd>& Rho, const Leaf& leaf, const Tree& tree);
template Matrix<d>
TreeIO::LeafDensity(const TensorTree<d>& Psi, const SparseMatrixTree<d>& Rho, const Leaf& leaf, const Tree& tree);


