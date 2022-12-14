//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef INTEGRATORINTERFACE_H
#define INTEGRATORINTERFACE_H
#include "TreeClasses/HamiltonianRepresentation.h"

class IntegratorInterface {
public:
	IntegratorInterface(const Hamiltonian& H, const Tree& tree, complex<double> phase = 1.)
		: hRep_(H, tree), H_(H), tree_(tree), phase_(phase), dPsi_(tree) {}

	~IntegratorInterface() = default;

	Wavefunction Derivative(double t, const Wavefunction& Psi) {
		::Derivative(dPsi_, hRep_, t, Psi, H_, tree_, phase_);
		return dPsi_;
	}
private:
	const Hamiltonian& H_;
	const Tree& tree_;
	complex<double> phase_;
	HamiltonianRepresentation hRep_;
	Wavefunction dPsi_;

};

#endif //INTEGRATORINTERFACE_H
