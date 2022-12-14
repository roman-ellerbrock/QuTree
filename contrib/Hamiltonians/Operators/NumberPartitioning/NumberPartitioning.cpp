//
// Created by Roman Ellerbrock on 2019-07-12.
//

#include "NumberPartitioning.h"

SOP SumNumbers(const vector<double>& n, const mctdhBasis& basis) {
	/*!
	 * \brief This creates the SOP S = \sum_i (n(i) * s_z(i)).
	 *
	 * This SOP can be squared to obtain the required hamiltonian for the
	 * number partitioning problem formulated as an quantum Ising model.
	 * Via convention, the last qubit (f+1) is set to the |1> and excluded
	 * from the simulation. Therefore, there should be one qubit less than
	 * numbers in the set n.
	 */
	size_t f = basis.nLeaves();
	assert(n.size() == f + 1);

	function<void(const PrimitiveBasis&, Tensorcd&, const Tensorcd&)> I = PauliMatrices::Identity;
	function<void(const PrimitiveBasis&, Tensorcd&, const Tensorcd&)> s = PauliMatrices::sigma_z;

	// S = \sum_i n(i)*s_z(i), where s_z(i) = pauli z-matrix, n(i) = i-th number of the set
	SOP S;
	for (size_t k = 0; k < f; ++k) {
		MultiParticleOperator M(s, k);
		S.push_back(M, n[k]);
	}
	// This sets the last number to group |1> manually (by choosing positive coefficient)
	MultiParticleOperator M(I, 0);
	S.push_back(M, n[f]);
	return S;
}

void NumberPartitioning::SpecialInitialize(const mctdhBasis& basis) {
	/*!
	 * \brief This creates the SOP that corresponds to the number-partitioning problem.
	 *
	 * This routine initializes the Hamiltonian for the number partitioning problem.
	 * The set of numbers are described by the vector<double> n. A primitive initialization
	 * for the numbers is implemented. The number of these numbers is adjusted to the number
	 * of qubits. The last number is asigned to the |1> state by freedom of choice and
	 * removed from the simulation.
	 */
	size_t f = basis.nLeaves();
	size_t fn = f + 1;
	// Create and print some numbers for the number partitioning problem
	// {1, 2, ..., N}
	vector<double> n;
	for (size_t k = 0; k < fn; ++k) {
		n.emplace_back(k + 1);
	}
	cout << "The numbers that are partitioned are:\t";
	for (auto x : n)
		cout << x << "\t";
	cout << endl;

	// H = (S)^2
	SOP S = SumNumbers(n, basis);
	S = S * S;

	for (size_t k = 0; k < S.size(); ++k) {
		push_back(S[k], S.Coeff(k));
	}

	cout << "Number of parts in Hamiltonian: " << size() << endl;
}
