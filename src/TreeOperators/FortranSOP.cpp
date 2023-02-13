#include "TreeOperators/FortranSOP.h"
#include "Core/Matrix.h"

extern "C"
{
	void DDXCPP(int* mode, double hPsi[], double Psi[], int* dim, double matrix[], double trafo[])
	{
		if (trafo != NULL) { return; }
		ApplyFortranDDX((complex<double>*) hPsi, (complex<double>*) Psi, (const LeafInterface&)*matrix, *dim);
	}

	void KINCPP(int* mode, double hPsi[], double Psi[], int* dim, double matrix[], double trafo[])
	{
		if (trafo != NULL) { return; }
		ApplyFortranKIN((complex<double>*) hPsi, (complex<double>*) Psi, (const LeafInterface&)*matrix, *dim);
	}

	void fillpara_(int* f, double* para);
}

void ApplyFortranDDX(complex<double> hPsi[], complex<double> Psi[],
	const LeafInterface& grid, int dim) {

	// Create a Tensordim
	TensorShape tdim({(size_t) dim, 1});

	// Create Tensors
	Tensorcd hPsivec(tdim);
	Tensorcd Psivec(tdim);
	for (int i = 0; i < dim; i++)
		Psivec(i) = Psi[i];

	// Apply p
	grid.applyP(hPsivec, Psivec);

	// p=-i*h*ddx -> ddx=i*p/h
	complex<double> imag(0., 1.);
	for (int i = 0; i < dim; i++)
		hPsi[i] = imag*hPsivec(i);
}

void ApplyFortranKIN(complex<double> hPsi[], complex<double> Psi[],
	const LeafInterface& grid, int dim) {
	// Create a Tensordim
	TensorShape tdim({(size_t) dim, 1});

	// Create Tensors
	Tensorcd hPsivec(tdim);
	Tensorcd Psivec(tdim);
	for (int i = 0; i < dim; i++)
		Psivec(i) = Psi[i];

	// Apply kin
	grid.applyKin(hPsivec, Psivec);

	// Copy back
	for (int i = 0; i < dim; i++)
		hPsi[i] = hPsivec(i);
}

FortranSOP::FortranSOP(const Tree& tree) {
	SpecialInitialize(tree);
}

void FortranSOP::callHinit(Vectorcd& coeffs, Matrix<int>& diag)
{
	// This routine should be overwritten for every inherited class
	// Call the Hinit function of this system

	cout << "This is only a template. Overwrite this function for a implementation." << endl;
	assert(0);

}

void FortranSOP::SpecialInitialize(const Tree& tree) {
//	hamiltonian.clear();
//	coeff.clear();
	size_t npart_coeff_start = coeff_.size();

	// Save the number of parts and modes in this hamiltonian
	InitOperator();

	// Initialize Fortran SysPar
	InitializeFortranPara(tree);

	// Call Hinit to get coefficients and diagonal terms in the hamiltonian
	Vectorcd fortrancoeffs(nparts + 1);
	Matrix<int> diag(nmodes, nparts);

	// Call Hinit with the new coefficients and diagonal-matrix
	callHinit(fortrancoeffs, diag);

	for (int i = 0; i < nparts + 1; i++) {
		if (fortrancoeffs(i) == 0.) {
			nparts = i;
			break;
		}
	}
	// cout << "nparts = " << nparts << endl;

	// cout << "Building Fortran Hamiltonian in C++..." << endl;
	// Push_back each Single-Particle Operator
	for (int i = 0; i < nparts; i++) {
		MLOcd M;
		for (int k = 0; k < nmodes; k++) {
			if (diag(k, i) != 1) {
				FortranOperator* F = new FortranOperator();
				const Leaf& phys = tree.getLeaf(k);
				int dim = phys.dim();
				int mode = phys.mode() + 1;
				F->Initialize(i + 1, mode, dim, SystemH);
				shared_ptr<LeafOperatorcd> F_shared(F);
				M.push_back(F_shared, k);
	//			cout << "part = " << i << " mode = " << k << " diag = " << diag(k, i) << "\n";
			}
			else {
	//			cout << "part = " << i << " mode = " << k << " diag = " << diag(k, i) << endl;
			}
		}
		push_back(M, fortrancoeffs(i));
	}

	// Assure that the sizes of all vectors fit together
	assert(size() == (nparts + npart_coeff_start) );

	{
		MLOcd M;
		for (int k = 0; k < nmodes; k++)
		{
			FortranOperator* F = new FortranOperator;
			const Leaf& phys = tree.getLeaf(k);
			int dim = phys.dim();
			F->Initialize(-1, k + 1, dim, SystemH);
			shared_ptr<LeafOperatorcd> F_shared(F);
			M.push_back(F_shared, k);
		}
		DividingSurface = M;
	}
	// Set Flux Operator
	{
		MLOcd M;
		for (int k = 0; k < nmodes; k++)
		{
			FortranOperator* F = new FortranOperator;
			const Leaf& phys = tree.getLeaf(k);
			int dim = phys.dim();
			F->Initialize(-2, k + 1, dim, SystemH);
			shared_ptr<LeafOperatorcd> F_shared(F);
			M.push_back(F_shared, k);
		}
		FluxOperator = M;
	}

//	cout << "Hamiltonian initialized." << endl;

}

void FortranSOP::InitOperator()
{
	cout << "Overwrite this function and set nparts and nmodes." << endl;
	assert(0);
}

Wavefunction FortranSOP::ApplyFlux(const Wavefunction& Psi, const Tree& basis) {
	// Apply the Flux operator F = [h, H]
	return FluxOperator.apply(Psi, basis);
}

Wavefunction FortranSOP::ApplyDividingSurface(const Wavefunction& Psi, const Tree& basis)
{
	// Apply the dividing surface oeprator h
	return DividingSurface.apply(Psi, basis);
}

void FortranSOP::InitializeFortranPara(const Tree& basis)
{
	// Copy Parameters into syspar in Fortran
	int f = basis.nLeaves();
	Matrixd para(4, f);
	for (int k = 0; k < f; k++)
	{
		const Leaf& phy = basis.getLeaf(k);
		const PhysPar par = phy.par();
		para(0, k) = par.omega();
		para(1, k) = par.r0();
		para(2, k) = par.wfr0();
		para(3, k) = par.wfOmega();
	}

	fillpara_(&f, &para(0, 0));
}
