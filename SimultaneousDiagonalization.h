#pragma once
#include "SingleParticleOperator.h"
#include "Vector.h"
#include "stdafx.h"
#include <random>
#include "JacobiRotationFramework.h"

using namespace JacobiRotationFramework;

class SimultaneousDiagonalization
{
public:
	SimultaneousDiagonalization();
	~SimultaneousDiagonalization();

	// Initialize Simultaneous Diagonalization
	void Initialization(vector<SPOcd>& A, double eps_);

	// Perform the Simultaneous Diagonalization
	void Calculate(vector<SPOcd>& A, SPOcd & trafo);
	
protected:
	// Perform a cycle of rotations over all matrices in A
	void JacobiRotations(vector<SPOcd>& A, SPOcd & trafo);

	// Measure off-Diagonality
	double MeasureOffDiagonals(const vector<SPOcd>& A);

	// Measure Diagonality
	double MeasureDiagonality(vector<SPOcd>& A);

	// Preconditioning of SD
	void InitialTransformation(vector<SPOcd>& A, SPOcd & trafo);

	int dim;
	int nmat;
	double eps;
};

