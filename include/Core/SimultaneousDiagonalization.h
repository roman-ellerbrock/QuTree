#pragma once
#include "Core/FactorMatrix.h"
#include "Vector.h"
#include "stdafx.h"
#include <random>
#include "Core/JacobiRotationFramework.h"

using namespace JacobiRotationFramework;

/**
 * \class SimultaneousDiagonalization
 * \ingroup Core
 * \brief This class performs a simulatneous diagonalization.
 *
 * Attempts to diagonalize a set of, potentially not commuting, matrices.
 */

class SimultaneousDiagonalization {
public:
	SimultaneousDiagonalization() = default;
	~SimultaneousDiagonalization() = default;

	// Initialize Simultaneous Diagonalization
	void Initialization(vector<FactorMatrixcd>& A, double eps_);

	// Perform the Simultaneous Diagonalization
	void Calculate(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo);

protected:
	// Perform a cycle of rotations over all matrices in A
	void JacobiRotations(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo);

	// Measure off-Diagonality
	double MeasureOffDiagonals(const vector<FactorMatrixcd>& A);

	// Measure Diagonality
	double MeasureDiagonality(vector<FactorMatrixcd>& A);

	// Preconditioning of SD
	void InitialTransformation(vector<FactorMatrixcd>& A, FactorMatrixcd& trafo);

	int dim;
	int nmat;
	double eps;
};

