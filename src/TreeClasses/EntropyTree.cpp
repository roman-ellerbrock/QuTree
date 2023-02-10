//
// Created by Roman Ellerbrock on 9/4/20.
//

#include "TreeClasses/EntropyTree.h"
#include "TreeClasses/TensorTreeFunctions.h"
#include "TreeClasses/SpectralDecompositionTree.h"

void EntropyTree::initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		attributes_.push_back(0.);
	}
}

void EntropyTree::calculate(const TensorTreecd& Psi, const Tree& tree) {
	initialize(tree);
	auto rho = TreeFunctions::contraction(Psi, tree, true);
	auto X = SpectralDecompositionTreecd(rho, tree);
	for (const Node& node : tree) {
		const Vectord& occ = X[node].second;
		double S = 0.;
		for (size_t i = 0; i < occ.dim(); ++i) {
			if (occ(i) > 0.) {
				S -= occ(i) * log(occ(i));
			}
		}
		operator[](node) = S;
	}
}

void EntropyTree::perplexity(const TensorTreecd& Psi, const Tree& tree,
	size_t power) {
	calculate(Psi, tree);

	for (const Node& node : tree) {
		double s = operator[](node);
		operator[](node) = pow(exp(s), power);
	}

}

void EntropyTree::print(const Tree& tree) {
	for (const Node& node : tree) {
		node.info();
		cout << operator[](node) << endl;
	}
}

double EntropyTree::metric(const Tree& tree) const {
	double acc = 0.;
	for (const Node& node : tree) {
		double t = operator[](node);
		if (t > acc) { acc = t; }
	}
	return acc;
}



