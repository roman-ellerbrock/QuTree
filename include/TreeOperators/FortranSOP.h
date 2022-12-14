#pragma once
#include "TreeOperators/Hamiltonian.h"
#include "Core/Vector.h"
#include "TreeOperators/FortranOperator.h"
#include "Core/Matrix.h"

// Wrapper functions for the primitive Fortran-operators. 
void ApplyFortranDDX(complex<double> hPsi[], complex<double> Psi[], const LeafInterface& phys, int dim);
void ApplyFortranKIN(complex<double> hPsi[], complex<double> Psi[], const LeafInterface& phys, int dim);

typedef TensorTreecd Wavefunction;

class FortranSOP :
	public SOPcd
{
	// Wrapper for Fortran-Code. Class is used to link C++ and Fortran-Code.
	// If you want a new System with a Fortran-Operator, create a new class and inherit from this class.
	// Then overwrite the functions callHinit and InitOperator. See "CH3_quasiexact"-class for an example.
	// Please ensure to keep code portable! :-)
public:
	FortranSOP(const Tree& basis);
	FortranSOP() = default;
	~FortranSOP() = default;

	// Overwritten member-function of Hamiltonian. Apply Fluxoperator
	Wavefunction ApplyFlux(const Wavefunction & Psi, const Tree& tree);

	// Overwritten member-function of Hamiltonian. Apply Dividing surface
	Wavefunction ApplyDividingSurface(const Wavefunction & Psi, const Tree& tree);

	// Initialize Fortran SysPar block
	void InitializeFortranPara(const Tree& tree);

protected:
	void SpecialInitialize(const Tree& tree);

	// Number of parts in Fortran-Hamiltonian
	int nparts;
	// Number of physical coordinates in the operator
	int nmodes;

	vector<Matrixd> matrix, trafo;

	MLOcd FluxOperator;
	MLOcd DividingSurface;

	// Fortran-Hamiltonian
	FortranSystemH SystemH;

private:
	// Template function to call Hinit
	virtual void callHinit(Vectorcd& coeffs, Matrix<int>& diag);

	// Template function to set nparts, nmodes and SystemH
	virtual void InitOperator();
};

