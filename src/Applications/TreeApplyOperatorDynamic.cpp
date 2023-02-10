//
// Created by Roman Ellerbrock on 2/20/20.
//

#include "Applications/TreeApplyOperatorDynamic.h"
#include "TreeClasses/TreeIO.h"
#include "Util/RandomProjector.h"
#include "TreeClasses/EntropyTree.h"

namespace TreeFunctions {

	void applyOperator(TensorTreecd& Psi, Tree& tree,
		const SOPVectorcd& sop, const SOPVectorcd& sop_adj, double eps,
		size_t max_spf, size_t n_plus) {
		assert(sop.size() == sop_adj.size());

		microseconds time(0);
		high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

//		size_t xreg = 5;
//		string name("perp.l" + to_string(xreg) + ".a2.txt");
//		cout << "will write on file:" << name << endl;
//		getchar();
//		ofstream os(name, ios::app);

//		size_t chunk = sop.size() / xreg;

		for (size_t i = 0; i < sop.size(); ++i) {
			applyOperator(Psi, tree, sop[i], sop_adj[i], eps, max_spf, n_plus);

/*			if (sop.size() > xreg && ((i + 1) % chunk) == 0) {
				EntropyTree S;
				S.Perplexity(Psi, tree);
//				S.print(tree);
				auto p = S.metric(tree);
				os << p << endl;
			}*/

			t2 = chrono::high_resolution_clock::now();;
			TreeIO::statusTime(i, sop.size(), 1, 20, t1, t2, time);
			t1 = t2;
		}

		microseconds global_time = chrono::duration_cast<chrono::microseconds>(t2 - t1);
		cout << "Total time at end: " << to_string(global_time.count() / 1000000.) << "s\n";
		cout << "Time per operation: " << to_string(global_time.count() / (1000. * sop.size()) + 1) << "ms\n.";
		cout.flush();
	}

	void applyOperator(TensorTreecd& Psi, Tree& tree,
		const SOPcd& sop, const SOPcd& sop_adj,
		double eps, size_t max_spf, size_t n_plus) {
		// This routine performs the
		// || H|Psi>-|Psi'> || = 0 optimization
		// tr(H|Psi><Psi|H)/k on every bottom node
		// Build H2 Holes/Matrices
		cerr << "add rho to this.\n";
		shared_ptr<SparseTree> active = make_shared<SparseTree>(SparseTree(sop, tree));
//		vector<SparseMatrixTreecd> Holes = BuildHHoles(Psi, sop, sop_adj, tree, active);
		vector<SparseMatrixTreecd> Holes;
		vector<SparseMatrixTreecd> Hmats(sop.size(), SparseMatrixTreecd(active, tree));

		for (const Node *node_p : *active) {
			const Node& node = *node_p;
			// Merge the SPFs of the Psis
			// First represent them in the tree of the optimized new SPF set (for upper layers)
			if (!node.isToplayer()) {
				Node& parent = tree.getNode(node.parent().address());
				Node& node_x = tree.getNode(node.address());
				Psi[node] = applyLowerDynamic(Psi[node], sop, Hmats,
					Holes, node_x, parent, Psi[parent], eps, max_spf, n_plus);
			} else {
				Psi[node] = applyTop(Psi[node], Hmats, sop, node);
			}
		}
	}

	Tensorcd applyLowerDynamic(Tensorcd Phi, const SOPcd& sop,
		vector<SparseMatrixTreecd>& Hmats, const vector<SparseMatrixTreecd>& Holes,
		Node& node, Node& parent, Tensorcd& upPhi, double eps, size_t max_spf, size_t n_plus) {
		vector<Tensorcd> Phis;
		for (size_t n = 0; n < sop.size(); ++n) {
			Phis.push_back(apply(Hmats[n], Phi, sop(n), node)); // @TODO: make Hmats constant here
		}

		// Calculate Eigenvectors/Eigenvalues of the "Delta"-operator and resulting new SPFs
		Matrixcd delta = buildDeltaOperator(Phis, sop, Holes, node);

		Matrixcd trafo(delta);
		Vectord eigenvalues(delta.dim1());
		delta.cDiag(trafo, eigenvalues);

		/*
		size_t p = 3;
		size_t rank = min(node.shape().lastDimension() * 3, node.shape().lastBefore());
		size_t thresh = node.shape().lastBefore();
		if (true) {
//		if (rank > node.shape().lastBefore() / 2) {
			delta.cDiag(trafo, eigenvalues);
		} else {
			mt19937 gen(time(NULL));
			auto x = Random::DiagonalizeRandom<complex<double>, Matrixcd>(delta, rank, p, gen);
			trafo = x.first;
			eigenvalues = x.second;
//			eigenvalues.print();
		}
		 */

		// Get number of occupied SPop
		// ts
		size_t n_occ = nOccupied_local(eigenvalues, eps);
		n_occ = max((size_t) 1, n_occ);
		double eps_warn = eps * 10;
		assert(n_occ > 0);
		n_occ += (size_t) n_plus;
		n_occ = min(n_occ, max_spf);
		n_occ = min(n_occ, node.shape().lastBefore());
		if (n_occ == 0) {
			node.info();
			cout << n_occ << endl;
			cout << node.shape().lastBefore() << endl;
			eigenvalues.print();
			exit(1);
		}

		/// Check whether max-cutting number of SPF violates hard accuracy limit
		if (n_occ == max_spf && n_occ != node.shape().lastBefore()) {
			if (eigenvalues(n_occ - 1) > eps_warn) {
				cerr << "WARNING: Critical accuracy cut-off due to maximum number of SPFs ("
					 << eigenvalues(n_occ - 1) << ").\n";
			}
		}

		// Adjust dim of Phi (no need to set SPFs, I think)
		Phi = Phi.adjustStateDim(n_occ);
		upPhi = upPhi.adjustActiveDim(n_occ, node.childIdx());
		TensorShape& tdim = node.shape();
		tdim.setDimension(n_occ, tdim.lastIdx());
		TensorShape& tdim_p = parent.shape();
		tdim_p.setDimension(n_occ, node.childIdx());
		assert(node.shape().lastDimension() == n_occ);
		assert(parent.shape()[node.childIdx()] == n_occ);

		// Create the new Tensor from the eigenvectors
		Tensorcd mPhi = occupyFromMatrix(trafo, node.shape());

		// Calculate and store overlaps of new SPFs with SPFs of Phis
		WorkMemorycd mem(node);
		for (size_t n = 0; n < sop.size(); ++n) {
			representLayer(Hmats[n], mPhi, Phi, sop(n), node, &mem);
		}

		return mPhi;
	}

	size_t nOccupied_local(const Vectord& ev, double rho) {
		size_t n_occ = ev.dim();
		for (int i = ev.dim() - 1; i >= 0; --i) {
			if (ev(i) < rho) {
				return ev.dim() - 1 - i;
			}
		}
		return n_occ;
	}
}
