#pragma once
#include "Tensor/Tensor.h"
#include "stdafx.h"
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

namespace simultaneousDiagonalization {

	// Perform the Simultaneous Diagonalization
	void calculate(vector<Matrixcd>& A, Matrixcd& trafo, double eps = 1e-10);

	// Perform a cycle of rotations over all matrices in A
	void jacobiRotations(vector<Matrixcd>& A, Matrixcd& trafo);

	// Measure off-Diagonality
	double measureOffDiagonals(const vector<Matrixcd>& A);

	// Measure Diagonality
	double measureDiagonality(vector<Matrixcd>& A);

	// Preconditioning of SD
	void initialTransformation(vector<Matrixcd>& A, Matrixcd& trafo);

};

