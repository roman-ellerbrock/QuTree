//
// Created by Roman Ellerbrock on 12/16/20.
//
#include "Applications/ApplyFirstOrder.h"
#include "Util/RandomProjector.h"
#include "TreeClasses/EntropyTree.h"
#include "TreeClasses/SparseTree.h"
#include "Applications/TreeApplyOperator.h"
#include "TreeClasses/TreeIO.h"

using namespace chrono;

namespace TreeFunctions {

	void applyOperatorSCF(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPVectorcd& sop, const Tree& tree) {

		microseconds time(0);
		high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

		ofstream os("fidelity.dat");
		os << "# p = 2; 2^6 layout with chi = 16" << endl;
		double f_av, F = 1., logF = 0.;
		for (size_t i = 0; i < sop.size(); ++i) {
			double f = applyOperatorIteration(Psi, rho, sop[i], tree);

			/// Evaluate fidelity & corresponding statistics
			F *= f;
			logF += log(f);
			f_av = exp(logF / (double) (i + 1));

			t2 = chrono::high_resolution_clock::now();;
			TreeIO::statusTime(i, sop.size(), 1, 20, t1, t2, time);
			cout << "\tf = " << f << " | F = " << F << " | f_geoav = " << f_av << flush;
			os << i << "\t" << f << "\t" << F << "\t" << f_av << endl;
			t1 = t2;
		}
		cout << endl;
	}

/// apply a SOPcd operator onto a wavefunction (the sop operator is needed in its adjoint form as well)
	double applyOperatorIteration(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPcd& sop, const Tree& tree) {
		/**
		 * This routine performs the
		 * || |Psi'> - H |Psi> || = 0 optimization.
		 * */

		if (sop.size() == 1) {
			sop[0].applyReference(Psi, tree);
			return 1.;
		}

		auto stree = make_shared<SparseTree>(sop, tree, false);

		/// Generate random tensors
		mt19937 gen(time(nullptr));
		TensorTreecd HPsi(gen, tree, false);
		for (const Node& node: tree) {
			if (!stree->isActive(node)) {
				HPsi[node] = Psi[node];
			}
		}

		size_t n_iter = 100;
		double f_last = 1.;
		double f = 0.;
		WorkMemorycd mem(tree);
		for (size_t i = 0; i < n_iter; ++i) {
			SparseMatrixTreescd Hmats = TreeFunctions::represent(sop, HPsi, Psi, stree, tree, &mem);
			SparseMatrixTreescd HHoles = TreeFunctions::contraction(HPsi, Psi, Hmats, rho, stree, tree);
			applyOperatorIteration(HPsi, Psi, rho, Hmats, HHoles, sop, *stree, tree);
			f_last = f;
//			f = fidelity(HPsi, Psi, Hmats, rho, sop, *stree, tree);
			f = fidelity(HPsi, Psi, sop, tree);
			if (((1.-f) < 1e-6) || (f - f_last) < 1e-10) { break; }
		}
		Psi = HPsi;
		TreeFunctions::contraction(rho, Psi, *stree, true);
		return f;

	}

	void applyOperatorIteration(TensorTreecd& HPsi, const TensorTreecd& Psi, MatrixTreecd& rho,
		SparseMatrixTreescd& Hmats, SparseMatrixTreescd& HHoles,
		const SOPcd& sop, const SparseTree& stree, const Tree& tree) {

		WorkMemorycd mem(tree);
		for (const Node* node_ptr : stree) {
			const Node& node = *node_ptr;
			bool top = node.isToplayer();
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				/// If parent node is not active, the equations simplify to toplayer equations
				if (!(stree.isActive(parent))) { top = true; }
			}
			/// @TODO: move 3 lines down and test
			HPsi[node] = applyOperatorLocal(Psi[node], rho, Hmats, HHoles, sop, node);
			if (!top) {
				/// Represent result in new basis
				for (size_t i = 0; i < sop.size(); ++i) {
					TreeFunctions::representLayer(Hmats[i], HPsi[node], Psi[node], sop(i), node, &mem);
				}
			} else {
				HPsi[node] = applyTop(Psi[node], Hmats, sop, node);
				gramSchmidt(HPsi[node]);
			}
		}
	}

	Tensorcd applyOperatorLocal(const Tensorcd& Phi, const MatrixTreecd& rho,
		const SparseMatrixTreescd& mats, const SparseMatrixTreescd& holes,
		const SOPcd& sop, const Node& node) {
		assert(mats.size() == holes.size());

		Tensorcd HPhi(Phi.shape());
		for (size_t i = 0; i < mats.size(); ++i) {
			Tensorcd mPhi = TreeFunctions::apply(mats[i], Phi, sop(i), node);
			const SparseMatrixTreecd& hole = holes[i];
			if (hole.isActive(node)) {
				mPhi = matrixTensor(hole[node], mPhi, node.parentIdx());
			} else if (!node.isToplayer()) {
				mPhi = matrixTensor(rho[node], mPhi, node.parentIdx());
			}
			HPhi += mPhi;
		}
//		GramSchmidt(HPhi);
		HPhi = qr(HPhi);

		return HPhi;
	}

	double fidelity(const TensorTreecd& Psi, const TensorTreecd& Psilast,
		const SparseMatrixTreescd& umat, const MatrixTreecd& rho,
		const SOPcd& u, const SparseTree& stree, const Tree& tree) {

		const Node& top = stree.node(stree.size() - 1);
		auto UPsi = applyTop(Psilast[top], umat, u, top);
		if (!top.isToplayer()) {
			UPsi = matrixTensor(rho[top], UPsi, top.parentIdx());
		}
		auto fmat = Psi[top].dotProduct(UPsi);

		complex<double> fam = fmat(0, 0);
		return pow(abs(fam), 2);
	}

}
