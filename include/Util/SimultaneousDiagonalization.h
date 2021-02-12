#pragma once
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "Core/stdafx.h"
#include <random>
#include "JacobiRotationFramework.h"

using namespace JacobiRotationFramework;

/**
 * \class SimultaneousDiagonalization
 * \ingroup Util
 * \brief This class performs a simulatneous diagonalization.
 *
 * Attempts to diagonalize a set of, potentially not commuting, matrices.
 */

class SimultaneousDiagonalization {
public:
	SimultaneousDiagonalization() = default;
	~SimultaneousDiagonalization() = default;

	// Initialize Simultaneous Diagonalization
	void initialization(vector<Matrixcd>& A, double eps);

	// Perform the Simultaneous Diagonalization
	void calculate(vector<Matrixcd>& A, Matrixcd& trafo);

protected:
	// Perform a cycle of rotations over all matrices in A
	void jacobiRotations(vector<Matrixcd>& A, Matrixcd& trafo);

	// Measure off-Diagonality
	double measureOffDiagonals(const vector<Matrixcd>& A);

	// Measure Diagonality
	double measureDiagonality(vector<Matrixcd>& A);

	// Preconditioning of SD
	void initialTransformation(vector<Matrixcd>& A, Matrixcd& trafo);

	int dim_;
	int nmat_;
	double eps_;
};

