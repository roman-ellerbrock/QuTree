//
// Created by Roman Ellerbrock on 4/10/20.
//

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H
#include "TreeClasses/TensorTree.h"
#include <random>
#include "TreeClasses/SparseMatrixTreeFunctions.h"

namespace FullRank {
	typedef Tensorcd Wavefunction;
}

namespace Measurements {

	typedef vector<size_t> Measurement;
	typedef map<Measurement, size_t> Sample;
	typedef vector<pair<vector<size_t>, size_t>> SampleVector;

	void Update(SparseMatrixTreecd& S, SparseMatrixTreecd& rho,
		const TensorTreecd& Psi,
		size_t next, size_t last, const Tree& tree);

	Measurement measurement(TensorTreecd& Psi, mt19937& gen,
		const vector<size_t>& targets, const Tree& tree);

	Vectord probabilityDistibution(const FullRank::Wavefunction& Psi);

	Sample sample(const TensorTreecd& Psi, mt19937& gen,
		size_t n_samples, const vector<size_t>& targets, const Tree& tree);

	Sample sample(const FullRank::Wavefunction& Psi, mt19937& gen,
		size_t number_samples);

	/// expected values of x, E_p(x)
	double expectX(const Vectord& p, const Vectord& x);

	/// expected values for moments of x, E_p(x^moment)
	double expectMomentX(const Vectord& p, Vectord x, size_t moment);

	///////////////////////////////////////////////////////////////
	/// Entropy-based measurements
	///////////////////////////////////////////////////////////////

	double crossEntropy(const FullRank::Wavefunction& p,
		const FullRank::Wavefunction& q);

	double crossEntropy(const Sample& X, const FullRank::Wavefunction& Psi,
		bool verbose);

	double kullbackLeiblerDivergence(const FullRank::Wavefunction& p,
		const FullRank::Wavefunction& q);

	double klDivergence(const Sample& X, const FullRank::Wavefunction& Psi,
		const function<double(const FullRank::Wavefunction& Psi)>& s, bool verbose);

	double crossEntropyDifference(const FullRank::Wavefunction& p,
		const FullRank::Wavefunction& q);

	double entropy(const FullRank::Wavefunction& Psi);

	double sPorterThomas(size_t n_qubits);

	void print(const Measurement& M, ostream& os);
	void print(const Sample& S, ostream& os = cout);

	ostream& operator<<(ostream& os, const Measurement& M);
}

#endif //MEASUREMENTS_H
