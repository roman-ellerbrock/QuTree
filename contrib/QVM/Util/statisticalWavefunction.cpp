//
// Created by Roman Ellerbrock on 2/19/21.
//

#include "statisticalWavefunction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/TreeIO.h"
#include "Circuits/GateOperators.h"

namespace statistical {

	Wavefunctions read(const string& name, const string& suffix) {
		Wavefunctions psi;
		if (!fileExists(name + suffix)) { return psi; }
		size_t i = 0;
		while (true) {
			if (!fileExists(name + to_string(i++) + suffix)) {
				return psi;
			}
			TensorTreecd npsi(name);
			psi.emplace_back(npsi);
		}
	}

	void Update(spMatrixTrees& S, spMatrixTrees& rho,
		const Wavefunctions& Psi,
		size_t next, size_t last, const Tree& tree) {

		for (size_t i = 0; i < Psi.size(); ++i) {
			Measurements::Update(S[i], rho[i], Psi[i], next, last, tree);
		}
	}

	Matrixcd leafDensity(const Wavefunctions& Psi, const spMatrixTrees& rho,
		const Leaf& leaf, const Tree& tree) {
		Matrixcd p(leaf.dim(), leaf.dim());
		for (size_t i = 0; i < Psi.size(); ++i) {
			p += TreeIO::leafDensity(Psi[i], rho[i], leaf, tree);
		}
		p /= p.trace();
		return p;
	}

	void projectPrimitive(Wavefunctions& Psis, size_t state, size_t coord, const Tree& tree) {
		/// Project on single-qubit state
		shared_ptr<LeafOperatorcd> P = make_shared<Circuits::PrimitiveProjector>(state);
		MLOcd Proj;
		Proj.push_back(P, coord);
		for (TensorTreecd& Psi : Psis) {
			Psi = Proj.apply(Psi, tree);
		}
	}

	size_t measureQubit(mt19937& gen, const Matrixcd& rho) {

		uniform_real_distribution<double> dist(0., 1.);
		double random_double = dist(gen);
		double acc = 0.;
		size_t i = 0;
		for (; i < rho.dim1(); ++i) {
			acc += real(rho(i, i));
			if (acc >= random_double) {
				break;
			}
		}
		if (i >= rho.dim1()) {
			cerr << "Error: measurement failed critically.\n";
			exit(1);
		}
		return i;
	}

	Measurement measurement(Wavefunctions& Psis, mt19937& gen,
		const vector<size_t>& targets, const Tree& tree) {

//		SparseMatrixTreecd overlap(targets, tree);
//		SparseMatrixTreecd holeOverlap(targets, tree);
		spMatrixTrees overlap;
		for (const auto& psi : Psis) { overlap.emplace_back(targets, tree); }
		spMatrixTrees holeOverlap;
		for (const auto& psi : Psis) { holeOverlap.emplace_back(targets, tree); }
		uniform_real_distribution<double> dist(0., 1.);
		Measurement measures;

		assert(!targets.empty());
		size_t idx = 0;
		for (size_t coord : targets) {
			const Leaf& leaf = tree.getLeaf(coord);
			const Node& node = (Node&) leaf.parent();

			/// calculate occupancy of single-qubit state (|0>, |1>)
			size_t last = targets[idx];
			size_t next = targets[idx];
			if (idx + 1 < targets.size()) { next = targets[idx + 1]; }
			Update(overlap, holeOverlap, Psis, next, last, tree);
			auto occupancy = leafDensity(Psis, holeOverlap, leaf, tree);

			size_t state = measureQubit(gen, occupancy);
			measures.emplace_back(state);

			/// Project on single-qubit state
			projectPrimitive(Psis, state, coord, tree);

			idx++;
		}
		return measures;
	}




}
