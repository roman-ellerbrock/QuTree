//
// Created by Roman Ellerbrock on 3/1/20.
//

#ifndef HAMILTONIANREPRESENTATION_H
#define HAMILTONIANREPRESENTATION_H
#include "TreeOperators/Hamiltonian.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Util/QMConstants.h"
#include "DVR/CDVR.h"
#include "DVR/MatrixTensorTreeFunctions.h"

#include "TreeOperators/ActiveCounter.h"

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation(const Hamiltonian& H, const Tree& tree,
		const Tree& cdvrtree, bool tail = true)
		: rho_(tree), rho_decomposition_(tree), hCorr_(tree), hConCorr_(tree),
		  rho_inverse_(tree), cdvr_(tree), mem_(tree), corrections_(!tail) {

		hMats_.clear();
		hContractions_.clear();
		size_t l = 0;
		for (const auto& M: H) {
			hMats_.emplace_back(SparseMatrixTreecd(M, tree, tail));
			hContractions_.emplace_back(SparseMatrixTreecd(M, tree, tail));
		}

		nActives_.initialize(hMats_, hContractions_, H, tree);
	}

	HamiltonianRepresentation(const Hamiltonian& H, const Tree& tree)
		: HamiltonianRepresentation(H, tree, tree) {
	}

	void initializeDense(const Hamiltonian& H, const Tree& tree) {
		/**
		 * \brief initialize with h-mat and h-hole matrices at every node.
		 */

		/// Create a vector of all modes
		vector<size_t> modes = {};
		for (size_t i = 0; i < tree.nLeaves(); ++i) {
			modes.push_back(tree.getLeaf(i).mode());
		}

		/// Create the hMats & hContraction with all modes
		hMats_.clear();
		hContractions_.clear();
		for (const auto& M: H) {
			hMats_.emplace_back(SparseMatrixTreecd(modes, tree, true));
			hContractions_.emplace_back(SparseMatrixTreecd(modes, tree));
		}

	}

	~HamiltonianRepresentation() = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree,
		double time) {

		/// Calculate density matrix tree
		TreeFunctions::contraction(rho_, Psi, tree, true);

		/// Density matrix tree decomposition
		rho_decomposition_.calculate(rho_, tree);

		/// Calculate inverse density matrix tree
		rho_inverse_ = rho_decomposition_.invert(tree, 1e-10);
//		rho_inverse_ = inverse(rho_, tree);

		/// Calculate h-matrix trees
		TreeFunctions::represent(hMats_, H, Psi, Psi, tree);

		/// Calculate h-matrix tree contractions
		TreeFunctions::contraction(hContractions_, hMats_, Psi, Psi, tree);

		/// build correction matrices
		if (corrections_) { buildCorrection(H, Psi, tree); }

		/// Calculate CDVR
//		ofstream os("points.dat");
		if (H.hasV) { cdvr_.Update(Psi, H.V_, tree, 0, false); }
	}

	void buildUp(const Hamiltonian& H, const Wavefunction& Psi,
		const Node& node, double time) {
		for (size_t l = 0; l < hMats_.size(); ++l) {
			TreeFunctions::representLayer(hMats_[l], Psi[node], Psi[node], H[l], node, &mem_);
		}

		/// build corrections
		if (corrections_) {
			buildUpCorrection(H, Psi, node);
		}

		if (H.hasV) {
			cerr << "Add local building of CDVR matrices to HamiltonianRepresentation.\n";
			exit(1);
		}
	}

	void buildDown(const Hamiltonian& H, const Wavefunction& Psi,
		const Node& node, double time) {
		const MatrixTreecd* null = nullptr;

		if (!node.isToplayer()) {
			for (size_t l = 0; l < hContractions_.size(); ++l) {
				const SparseTree& stree = hContractions_[l].sparseTree();
				if (stree.isActive(node)) {
					TreeFunctions::contractionLayer(hContractions_[l], Psi, Psi, hMats_[l],
						null, stree, node, &mem_);
				}
			}
		}

		/// build corrections
		if (corrections_) {
			buildDownCorrection(H, Psi, node);
		}

		if (H.hasV) {
			cerr << "Add local building of CDVR matrices to HamiltonianRepresentation.\n";
			exit(1);
		}
	}

	void buildSCF(const Hamiltonian& H, const Wavefunction& Psi,
		const Tree& tree, double time) {

		/// Calculate h-matrix trees
//		TreeFunctions::represent(hMats_, H, Psi, Psi, tree);
		for (const Node& node : tree) {
			buildUp(H, Psi, node, time);
		}

		/// Calculate h-matrix tree contractions
//		TreeFunctions::contraction(hContractions_, hMats_, Psi, Psi, tree);
		for (int n = tree.nNodes()- 1; n >= 0; --n) {
			const Node& node = tree.getNode(n);
			buildDown(H, Psi, node, time);
		}

		if (corrections_) { buildCorrection(H, Psi, tree); }
	}

	void buildSCF(const Hamiltonian& H, const Wavefunction& Psi,
		const Node& node, const Node& next, double time) {

		/// only works for sane node & next input!
		bool up = true;
		if (!node.isToplayer()) {
			/// node.parent == next
			if (node.parent().address() == next.address()) {
				up = true;
			} else {
				up = false;
			}
		} else {
			/// next.parent == next
			if (next.parent().address() == node.address()) {
				up = false;
			} else {
				up = true;
			}
		}

		if (up) {
			buildUp(H, Psi, node, time);
		} else {
			buildDown(H, Psi, next, time);
		}

	}

	void buildUpCorrection(const Hamiltonian& H, const Wavefunction& Psi, const Node& node);
	void buildDownCorrection(const Hamiltonian& H, const Wavefunction& Psi, const Node& node);
	void buildCorrection(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree);

	void print(const Tree& tree, ostream& os = cout) {
		os << "Rho:" << endl;
		rho_.print(tree, os);
		os << "Matrix representations:" << endl;
		for (const auto& mat: hMats_) {
			mat.print(os);
		}
		os << "Matrix Contractions:" << endl;
		for (const auto& con: hContractions_) {
			con.print(os);
		}
	}

	int counter_ = 0;

	MatrixTreecd rho_;
	SpectralDecompositionTreecd rho_decomposition_;
	MatrixTreecd rho_inverse_;
	SparseMatrixTreescd hMats_;
	SparseMatrixTreescd hContractions_;

	bool corrections_; /// decides whether corrections are turned on or off
	SparseMatrixTreecd hCorr_; /// correction term for hMats
	SparseMatrixTreecd hConCorr_; /// correction term for hCon

	ActiveCounter nActives_;

	WorkMemorycd mem_;

	CDVR cdvr_;
};

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep, const Node& node);

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree);

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase = 1.);

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase = 1.);

void symDerivative(MatrixTensorTree& dPsi, HamiltonianRepresentation& hRep,
	double time, const MatrixTensorTree& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase = 1.);

void output(const HamiltonianRepresentation& hrep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree);

#endif //HAMILTONIANREPRESENTATION_H
