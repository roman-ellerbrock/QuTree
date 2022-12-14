#include "TreeOperators/FortranOperator.h"

void FortranOperator::Initialize(int part_, int mode_, int dim_, FortranSystemH SystemH_) {
	// Set relevant fortran-variables:
	// Number of parts in the SOP-operator
	part = part_;
	// Number of physical coordinates
	mode = mode_;
	// Dimension of primitive grid
	dim = dim_;
	// Hamiltonian (Fortran-function)
	SystemH = SystemH_;
}

void FortranOperator::apply(const LeafInterface& grid, Tensorcd& hAcoeff, const Tensorcd& Acoeff)const {
	// Call the Fortran-Hamiltonian (using functional programming)
	// Intructions:
	// 1. Trafo and Matrix may NOT be used in the hamiltonian directly. So the routine SystemH (aka h)
	// is not allowed to operate with Matrix or Trafo. Usually the old mctdh-library used these routines
	// and there are now wrappers which take care of these operations in C++.
	// 2. Trafo is set to NULL. This fact is used in the library routines KIN, DDX, ... (See FortranHamiltonian.cpp).
	// The H+CH4-operator doesnt work anymore, if you set Trafo to another value.
	// Maybe also other operators! So dont change it, 
	// unless you know what you are doing.
	// 3. x, Psi, HPsi can be directly accessed in SystemH.
	Vectord x = grid.getX();
	for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
		SystemH((int*)&mode, (int*)&part, (double*)&hAcoeff(n*dim),
			(double*)&Acoeff(n*dim), (int*)&dim, (double*)&grid, NULL, &x(0));
	}

}
