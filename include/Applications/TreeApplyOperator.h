//
// Created by Roman Ellerbrock on 2/19/20.
//

#ifndef TREEAPPLYOPERATOR_H
#define TREEAPPLYOPERATOR_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SOPVector.h"

typedef TensorTreecd Wavefunction;

using namespace chrono;

/**
 * \brief The mlmcApply namespace contains functions that can be used to apply operators on a wavefunction.
 *
 * The namespace provides a method for applying SOPcd operators on wavefunctions
 * that scales quadratically with respect to the number of summands. Furthermore, some functions
 * are used by other classes that allow to apply SOPcd operators on wavefunctions.
 *
 * The present scheme for applying operators is described in Ref. [1]
 *
 * [1] R. Ellerbrock, T. J. Martínez, to be submitted (2019)
 */

struct Fidelity {
	/// simple tracker for fidelity

	void append(double f) {
		history_.emplace_back(f);
		total_ *= f;
		logTotal_ += log(f);
		++numberOperations_;
		avg_ = exp(logTotal_ / (double) (numberOperations_));
	}

	void calculate(const TensorTreecd& Psi, const TensorTreecd& Psilast,
		const SOPcd& u, const Tree& tree, WorkMemorycd* mem);

	void print(ostream& os = cout) const {
		os << "\tf = " << history_.back() << " | F = " << total_ << " | avg = " << avg_ << flush;
	}

	void write(const string& name = "fidelity.dat", const string& header = "") {
		/// Write history to file
		os_ = ofstream(name);
		auto hist = history_;
		history_.clear();
		avg_ = 1.;
		total_ = 1.;
		logTotal_ = 0.;
		numberOperations_ = 0;
		if (!header.empty()) { os_ << header << endl;}
		for (double f : hist) {
			append(f);
			os_ << history_.back() << "\t" << avg_ << "\t" << logTotal_ << endl;
		}
	}

	vector<double> history_;
	double avg_ = 1.; /// geometric average
	double total_ = 1.;
	double logTotal_ = 0.;
	size_t numberOperations_ = 0;
	ofstream os_ = ofstream("fidelity.dat");
};

namespace TreeFunctions {
	struct Diagonalizer;

	void applyOperator(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPVectorcd& sop, const SOPVectorcd& sop_adj,
		const Tree& tree, Fidelity& f, mt19937& gen);

	/// Routine for applying a SOPcd operator onto a wavefunction
	void applyOperator(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPcd& sop, const SOPcd& sop_adj, const Tree& tree, const Diagonalizer& par, WorkMemorycd* mem);

	/// Contract multiple wavefunctions in a single wavefunction
	TensorTreecd add(const vector<TensorTreecd>& Psis,
		const Tree& tree, double eps);

	double fidelity(const TensorTreecd& Psi, const TensorTreecd& Psilast,
		const SOPcd& u, const Tree& tree);

	/////////////////////////////////////////////////////////////////////
	// Helper functions
	/////////////////////////////////////////////////////////////////////

	/// Occupy a Tensor from a matrix
	Tensorcd occupyFromMatrix(const Matrixcd& ev, const TensorShape& tdim);

	/// Perform a state-index matrix-tensor operators for a set of Tensors/Matrices and accumulate
	Tensorcd multAdd(const vector<Tensorcd>& A, const vector<Matrixcd>& rho);

	/// Solve equations for applying an operator for a lower node
	Tensorcd applyLower(const Tensorcd& Phi, const SOPcd& sop,
		vector<SparseMatrixTreecd>& Hmats, const vector<SparseMatrixTreecd>& Holes,
		const Node& node, const Diagonalizer& par, WorkMemorycd* mem);

	/// Solve equations for applying an operator for the top node
	Tensorcd applyTop(const Tensorcd& Acoeff, const vector<SparseMatrixTreecd>& Hmats,
		const SOPcd& sop, const Node& node);

	/// Build hole-matrices used for applying scheme
	vector<SparseMatrixTreecd> buildHHoles(const TensorTreecd& Psi,
		const SOPcd& sop, const SOPcd& sop_adj, const MatrixTreecd& rho,
		const Tree& tree, shared_ptr<SparseTree>& active, WorkMemorycd* mem);

	/// Build residual operator from a set of Tensors, hole-matrices and the SOPcd operator.
	Matrixcd buildDeltaOperator(const vector<Tensorcd>& Phis, const SOPcd& sop,
		const vector<SparseMatrixTreecd>& Holes, const Node& node);

};

#endif //TREEAPPLYOPERATOR_H
