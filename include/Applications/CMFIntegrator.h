//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef CMFINTEGRATOR_H
#define CMFINTEGRATOR_H
#include "TreeClasses/HamiltonianRepresentation.h"
#include "TreeClasses/LayerInterface.h"
#include "TreeClasses/IntegratorVariables.h"
#include "Util/BS_integrator.h"

/**
 * \defgroup Integrator
 * \ingroup mctdh
 * \brief This module holds all classes related to mctdh-Integrators.
 *
 * */

/**
 * \class CMFIntegrator
 * \ingroup Integrator
 * \brief The constant mean-field integrator for solving the
 * mctdh-EOM.
 *
 * In the CMF-Integrator the mctdh-Matrices are assumed to be constant for
 * a given time to propagate a Wavefunction.
 *
 * */

typedef BS_integrator<LayerInterface&, Tensorcd, complex<double>> bs_integrator;

class CMFIntegrator {
public:
	CMFIntegrator(const Hamiltonian& H,
		const Tree& tree, const Tree& cdvrtree, complex<double> phase);

	CMFIntegrator(const Hamiltonian& H,
		const Tree& tree, complex<double> phase)
		: CMFIntegrator(H, tree, tree, phase) {}

	~CMFIntegrator() = default;

	void Integrate(IntegratorVariables& ivars, ostream& os = cout);

	double Error(const Wavefunction& Psi, const Wavefunction& Chi, const MatrixTreecd& rho, const Tree& tree) const;

	void Output(double time, const Wavefunction& Psi,
		const Wavefunction& Psistart, const Hamiltonian& H,
		const Tree& tree, ostream& os);

private:

	void CMFstep(Wavefunction& Psi, double time, double timeend, double accuracy_leaf, const Tree& tree);

	ofstream cdvr_file_;
	size_t cdvr_nr_;

	HamiltonianRepresentation matrices_;

	vector<bs_integrator> bs_integrators_;
	vector<double> dt_bs_;
	vector<LayerInterface> interfaces_;

	double tconst_;
	double max_increase_;
};


#endif //CMFINTEGRATOR_H
