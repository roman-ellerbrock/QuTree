//
// Created by Roman Ellerbrock on 2/19/21.
//

#ifndef STATISTICALWAVEFUNCTION_H
#define STATISTICALWAVEFUNCTION_H
#include "TreeClasses/TensorTreeFunctions.h"
#include "Measurements.h"
#include "Util/filenames.h"

namespace statistical {
	typedef vector<SparseMatrixTreecd> spMatrixTrees;

	class Wavefunctions: public vector<TensorTreecd> {
	public:
		Wavefunctions() = default;

		Wavefunctions(const string& name, const string& ending = ".dat") {
			ifstream is(name + ending);
			if (!is.good()) {
				cerr << "Error: cannot find statistical wavefunction.\n";
				exit(1);
			}
			push_back(TensorTreecd(is));
			size_t num = 1;
			while ((is = ifstream(name + "." + to_string(num++) + ending)) && is.good()) {
				push_back(TensorTreecd(is));
			}
		}

		Wavefunctions(const initializer_list<TensorTreecd>& ts) {
			for (const auto& x : ts) {
				emplace_back(x);
			}
		}

		~Wavefunctions() = default;
	};

	using Measurements::Measurement;

	Wavefunctions read(const string& name, const string& suffix);

	Measurement measurement(Wavefunctions& Psis, mt19937& gen,
		const vector<size_t>& targets, const Tree& tree);

	void Update(spMatrixTrees& S, spMatrixTrees& rho,
		const Wavefunctions& Psi,
		size_t next, size_t last, const Tree& tree);

	Matrixcd leafDensity(const Wavefunctions& Psi, const spMatrixTrees& rho,
		const Leaf& leaf, const Tree& tree);

	size_t measureQubit(mt19937& gen, const Matrixcd& rho);
}

#endif //STATISTICALWAVEFUNCTION_H
