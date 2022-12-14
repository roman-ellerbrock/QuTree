#include "CH3_quasiexact.h"

extern "C" {
	void hinitchhh_(double* coeff, int* diag, int* nmodes);
	void hchhh_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
}

CH3_quasiexact::CH3_quasiexact(const Tree& tree) {
	SpecialInitialize(tree);
}

void CH3_quasiexact::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag) {
	// Call Hinit here
	hinitchhh_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
}

void CH3_quasiexact::InitOperator() {
	// Set nparts and nmodes
	nparts = 58;
	nmodes = 6;

	// Set Hamiltonian
	SystemH = hchhh_;
}

