//
// Created by Roman Ellerbrock on 2/19/20.
//

#include "Applications/TreeApplyOperator.h"
#include "Core/Tensor_Extension.h"
#include "Util/RandomProjector_Implementation.h"
#include "TreeOperators/DeltaOperator.h"
#include "TreeClasses/TreeIO.h"

namespace TreeFunctions {

	struct Diagonalizer {
		Diagonalizer(mt19937& gen)
			: gen_(gen) {}

		size_t p = 10;
		mt19937& gen_;
	};

	using namespace TreeFunctions;

	double fidelity(const TensorTreecd& Psi, const TensorTreecd& Psilast,
		const SOPcd& u, const Tree& tree) {

		auto stree = make_shared<SparseTree>(u, tree);
		WorkMemorycd mem(tree);
		vector<SparseMatrixTreecd> umat = TreeFunctions::represent(u, Psi, Psilast, stree, tree, &mem);
		const Node& top = tree.topNode();

		auto UPsi = applyTop(Psilast[top], umat, u, top);
		auto fmat = Psi[top].dotProduct(UPsi);

		complex<double> fam = fmat(0, 0);
		return pow(abs(fam), 2);
	}

	void applyOperator(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPVectorcd& sop, const SOPVectorcd& sop_adj,
		const Tree& tree, Fidelity& f, mt19937& gen) {
		assert(sop.size() == sop_adj.size());
		microseconds time(0);
		high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

		Diagonalizer par(gen);

		TensorTreecd Psilast = Psi;
		WorkMemorycd mem(tree);
		for (int i = 0; i < sop.size(); ++i) {
			/// apply the operator
			applyOperator(Psi, rho, sop[i], sop_adj[i], tree, par, &mem);

			/// Evaluate fidelity & corresponding statistics
			f.calculate(Psi, Psilast, sop[i], tree, &mem);
			Psilast = Psi;

			/// Output
			t2 = chrono::high_resolution_clock::now();;
			TreeIO::statusTime(i+1, sop.size(), 1, 20, t1, t2, time);
			f.print(cout);
//			if (f.history_.back() > 1.01) { cout << "size:" << sop.size() << " " << i << endl; getchar(); }

			t1 = t2;
		}
		/// Final output
		microseconds global_time = chrono::duration_cast<chrono::microseconds>(t2 - t1);
//		cout << "\rTotal time at end: " << to_string(global_time.count() / 1000000.) << "s\t";
//		cout << "Time per operation: " << to_string(global_time.count() / ((double) (1000. * sop.size()) + 1.))
//			 << "ms\n";

//		f.print();
//		f.write();
		cout.flush();
	}

	/// apply a SOPcd operator onto a wavefunction (the sop operator is needed in its adjoint form as well)
	void applyOperator(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPcd& sop, const SOPcd& sop_adj, const Tree& tree, const Diagonalizer& par, WorkMemorycd* mem) {
		/**
		 * This routine performs the
		 * || (1 - P) H |Psi> || = 0 optimization.
		 * It is solved by finding the highest eigenvalue eigenvectors of
		 * tr(H|Psi><Psi|H)/k on every bottom node in a single bottom-up swipe.
		 * */

		/// Catch trivial case where sop has a single summand.
		if (sop.size() == 1) {
			sop[0].applyReference(Psi, tree);
//			TreeFunctions::contraction(rho, Psi, tree, true);
			return;
		}

		/** For high-dimensional systems and operators that only act on some nodes, it is important
		 * to only touch nodes that are active. First we build a SparseTree which allows to
		 * iteratre over the active nodes.
		 * */
		shared_ptr<SparseTree> stree = make_shared<SparseTree>(SparseTree(sop, tree, false));

		/// Build matrix representations of the operators in the SOPcd operator.
		vector<SparseMatrixTreecd> Holes = buildHHoles(Psi, sop, sop_adj, rho, tree, stree, mem);
		vector<SparseMatrixTreecd> Hmats(sop.size(), SparseMatrixTreecd(stree, tree));

		/// Swipe over the stree nodes and solve for the optimal SPFs. Psi will contain the result.
		for (const Node *node_p : *stree) {
			const Node& node = *node_p;
			// Merge the SPFs of the Psis
			// First represent them in the tree of the optimized new SPF set (for upper layers)
			bool top = node.isToplayer();
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				/// If parent node is not active, the equations simplify to toplayer equations
				if (!(stree->isActive(parent))) { top = true; }
			}
			if (top) {
				Psi[node] = applyTop(Psi[node], Hmats, sop, node);
				gramSchmidt(Psi[node]);
			} else {
				if (node.shape().lastBefore() == node.shape().lastDimension()) {
					for (size_t n = 0; n < sop.size(); ++n) {
						representLayer(Hmats[n], Psi[node], Psi[node], sop(n), node, mem);
					}
					continue;
				}
				Psi[node] = applyLower(Psi[node], sop, Hmats, Holes, node, par, mem);
			}
		}
		TreeFunctions::contraction(rho, Psi, *stree, true);
	}

	Tensorcd applyTop(const Tensorcd& Acoeff, const vector<SparseMatrixTreecd>& Hmats,
		const SOPcd& sop, const Node& node) {
		/// apply the parts of the SOPcd operator and accumulate the result.
		// @TODO: do accumulation in-place, without storing Tensors
		Tensorcd nAcoeff(Acoeff.shape());
		for (size_t n = 0; n < sop.size(); ++n) {
			nAcoeff += apply(Hmats[n], Acoeff, sop(n), node);
		}

		/// Re-orthonormalize the top layer coefficients
//		GramSchmidt(nAcoeff);
		return nAcoeff;
	}

	Tensorcd ApplyDelta(const Tensorcd& R, const vector<Tensorcd>& hPhi, const vector<SparseMatrixTreecd>& Holes,
		size_t npart, const Node& node) {
		Tensorcd Delta_R(R.shape());
		for (size_t l = 0; l < npart; ++l) {
			Matrixcd ml = hPhi[l].dotProduct(R);
			for (size_t o = 0; o < npart; ++o) {
				const SparseMatrixTreecd& hole = Holes[npart * l + o];
				Matrixcd Mlo = hole[node] * ml;
				Delta_R += matrixTensor(Mlo, hPhi[o], node.parentIdx());
			}
		}
		return Delta_R;
	}

	Tensorcd applyLower(const Tensorcd& Phi, const SOPcd& sop,
		vector<SparseMatrixTreecd>& Hmats, const vector<SparseMatrixTreecd>& Holes,
		const Node& node, const Diagonalizer& par, WorkMemorycd* mem) {
		/// For each part "h" in the sop operator, build h|Phi> and store the result in a vector

		/// Build the residual operator "Delta" (Eq. 1.10 in Ref. [1]) and obtain the eigenvectors
		size_t rank = node.shape().lastDimension();

		DeltaOperator deltaOp(Phi, Hmats, Holes, sop, node);
		DeltaMemory mem2(Phi, Hmats, sop, node);
		auto x = Random::diagonalizeRandom<complex<double>, DeltaOperator, DeltaMemory>(deltaOp,
			rank, par.p, par.gen_, &mem2);
		auto trafo = x.first;

		/// Create the new Tensor from the eigenvectors
		Tensorcd mPhi = occupyFromMatrix(trafo, node.shape());

		/// calculate and store overlaps of new SPFs with SPFs of Phis
		for (size_t n = 0; n < sop.size(); ++n) {
			representLayer(Hmats[n], mPhi, Phi, sop(n), node, mem);
		}
		return mPhi;
	}

	vector<SparseMatrixTreecd> buildHHoles(const TensorTreecd& Psi,
		const SOPcd& sop, const SOPcd& sop_adj, const MatrixTreecd& rho,
		const Tree& tree, shared_ptr<SparseTree>& active, WorkMemorycd* mem) {
		size_t size = sop.size();
		vector<SparseMatrixTreecd> Holes(size * size, SparseMatrixTreecd(active, tree));
		SparseMatrixTreecd hmat(active, tree);

		for (size_t l = 0; l < size; ++l) {
			for (size_t m = 0; m < size; ++m) {
				size_t idx = size * l + m;
				MLOcd MN = sop_adj(l) * sop(m);
				represent(hmat, MN, Psi, Psi, tree, mem);
//				contraction(Holes[idx], Psi, Psi, hmat, rho, tree);
				contraction(Holes[idx], Psi, Psi, hmat, &rho, *active, tree, mem);
			}
		}
		return Holes;
	}

	Matrixcd buildDeltaOperator(const vector<Tensorcd>& Phis, const SOPcd& sop,
		const vector<SparseMatrixTreecd>& Holes, const Node& node) {
		size_t dimpart = node.shape().lastBefore();
		Matrixcd delta(dimpart, dimpart);
		size_t size = sop.size();
		for (size_t l = 0; l < size; ++l) {
			vector<Matrixcd> holes;
			for (size_t m = 0; m < size; ++m) {
				const SparseMatrixTreecd& hole = Holes[size * l + m];
				holes.push_back(hole[node]);
			}
			Tensorcd B = multAdd(Phis, holes);
			Tensor_Extension::OuterProductAdd(delta, B, Phis[l]);
		}
		return delta;
	}

	Tensorcd multAdd(const vector<Tensorcd>& A, const vector<Matrixcd>& rho) {
		assert(A.size() == rho.size());
		const TensorShape& tdim = A[0].shape();
		Tensorcd B(tdim);
		for (size_t l = 0; l < A.size(); ++l) {
			B += multStateAB(rho[l], A[l]);
		}
		return B;
	}

	Tensorcd occupyFromMatrix(const Matrixcd& ev, const TensorShape& tdim) {
		// Occupy a tensor using the last n eigenvectors of a transformation matrix.
		// Note that there is no re-orthonormalization!
		assert(tdim.lastBefore() == ev.dim1());
		assert(tdim.lastDimension() <= ev.dim2());
		Tensorcd Phi(tdim);
		size_t last1 = tdim.lastDimension() - 1;
		size_t last2 = ev.dim2() - 1;
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			for (size_t i = 0; i < tdim.lastBefore(); ++i) {
				Phi(i, last1 - n) = ev(i, last2 - n);
			}
		}
		return Phi;
	}
}

void Fidelity::calculate(const TensorTreecd& Psi, const TensorTreecd& Psilast,
	const SOPcd& u, const Tree& tree, WorkMemorycd* mem) {

	auto stree = make_shared<SparseTree>(u, tree);
	vector<SparseMatrixTreecd> umat = TreeFunctions::represent(u, Psi, Psilast, stree, tree, mem);
	const Node& top = tree.topNode();

	auto UPsi = TreeFunctions::applyTop(Psilast[top], umat, u, top);
	auto fmat = Psi[top].dotProduct(UPsi);

	complex<double> fam = fmat(0, 0);
	append(pow(abs(fam), 2));
}

