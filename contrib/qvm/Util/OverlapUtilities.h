//
// Created by Roman Ellerbrock on 1/19/21.
//

#ifndef OVERLAPUTILITIES_H
#define OVERLAPUTILITIES_H
#include "yaml-cpp/yaml.h"
#include "TreeShape/Tree.h"
#include "Util/statisticalWavefunction.h"

namespace Utility {
	void statisticalWavefunctionOverlap(const Tree& tree);

	double calcXEB(const vector<double>& ps, size_t n_qubits);

	void xeb(const Tree& fr_tree, const Tree& xeb_tree, const string& file_fr, const string& file_b);

	void xeb_stat(ostream& os, mt19937& gen, const TensorTreecd& Psi,
		const statistical::Wavefunctions& Chis, const Tree& tree,
		size_t nsample);

	void wavefunctionOverlap(const Tree& tree, const vector<TensorTreecd>& Psi,
		const vector<TensorTreecd>& Chi);
}

#endif //OVERLAPUTILITIES_H
