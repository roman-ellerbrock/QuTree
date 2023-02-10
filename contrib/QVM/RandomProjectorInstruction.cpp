//
// Created by Roman Ellerbrock on 1/18/21.
//

#include "RandomProjectorInstruction.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "GateOperators.h"

void RandomProjectorInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho,
	QVMState& state, ostream& os) {
	/// Transform to canonical representation
	const Tree& tree = state.tree_;
	canonicalTransformation(Psi, tree);
	rho = TreeFunctions::contraction(Psi, tree, true);
	size_t q = q_;
	const Leaf& leaf = tree.getLeaf(q);
	const Node& node = (Node&) leaf.parent();
	Matrixcd dens = rho[node];
	dens /= dens.trace();

	/// Roll which state to project onto
	uniform_real_distribution<double> dist(0., 1.);
	double r = dist(state.gen_);
	double acc = 0.;
	size_t i = 0;
	for (; i < dens.dim1(); ++i) {
		acc += real(dens(i, i));
		if (acc >= r) {
			break;
		}
	}
	cout << endl << "i = " << i << endl;
	auto Psi_last = Psi;

	shared_ptr<LeafOperatorcd> P = make_shared<Circuits::PrimitiveProjector>(i);
	MLOcd Proj;
	Proj.push_back(P, q);
	Psi = Proj.apply(Psi, tree);
	qrOrthogonal(Psi, tree);

	auto fs = TreeFunctions::dotProduct(Psi, Psi_last, tree);
	double f = abs(fs[tree.topNode()](0,0));
	cout << "f = " << f << endl;
//	getchar();
	state.f_.append(f);

	orthonormal(Psi, tree);
	++state.pos_;
}

void RandomProjectorInstruction::apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) {
	cerr << "Not implemented yet.\n";
	exit(1);
}
