//
// Created by Roman Ellerbrock on 11/14/22.
//

#ifndef SCF_H
#define SCF_H
#include "TreeShape/Tree.h"
#include "TreeOperators/Hamiltonian.h"


vector<const Node*> scf_sweep(const Tree& tree);

struct SCF_parameters {
	size_t nIter{20};  /// how many iterations
	size_t nKrylov{5}; /// how large is the krylov space
	size_t nITP{0};    /// how many steps of imaginary time propagation
	double beta{1.};   /// how large should beta be?
	size_t output{0};/// output wavefunction every iteration
	double conversion{219474.6313705e0};

	TensorTreecd* psi{nullptr};
	const Hamiltonian* h{nullptr};
	const Tree* tree{nullptr};
};

struct KrylovSpace {
	KrylovSpace(const TensorShape& shape, size_t size) :
		space_(size, Tensorcd(shape)) {
	}

	KrylovSpace(vector<Tensorcd> space, SpectralDecompositioncd spectrum) :
		space_(move(space)), spectrum_(move(spectrum)) {}
	vector<Tensorcd> space_;
	SpectralDecompositioncd spectrum_;
};


void scf(SCF_parameters& par);

#endif //SCF_H
