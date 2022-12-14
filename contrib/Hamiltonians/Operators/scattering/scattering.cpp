#include <QMConstants.h>
#include "scattering.h"

namespace scatteringoperators {
	void Applyfx(Tensorcd& HPsi, const Tensorcd& Psi, const Vectord& x, function<double(double)> f) {
		assert(HPsi.Dim() == Psi.Dim());
		assert(HPsi.Dim().getdimpart() == x.Dim());
		for (int n = 0; n < Psi.Dim().getntensor(); n++)
			for (int i = 0; i < Psi.Dim().getdimpart(); i++)
				HPsi(i, n) = Psi(i, n) * f(x(i));
	}

	double eckhart(double x, double height, double xscale) {
//		return height*2./pow(exp(-(x-x0)/xscale)+exp((x-x0)/xscale),2);
//		return height*2./pow(exp(-(x-x0)/xscale)+exp((x-x0)/xscale),2);
		return height / pow(cosh(x / xscale), 2);
	}

	void CAP::Apply(Tensorcd& HPsi, const PrimitiveBasis& grid, const Tensorcd& Psi)const {
		assert(HPsi.Dim() == Psi.Dim());
		const Vectord& x = grid.GetX();
		for (int n = 0; n < Psi.Dim().getntensor(); n++) {
			for (int i = 0; i < Psi.Dim().getdimpart(); i++) {
				double y = abs(x(i));
				if (y > x0) {
					double quartic = pow((y - x0) / length, 4) * strength;
					quartic = min(quartic, strength);
					HPsi(i, n) = -Psi(i, n) * QM::im * quartic;
				} else {
					HPsi(i, n) = 0.;
				}
			}
		}
	}

	void V1(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi) {
		auto f = [](double x) { return eckhart(x, 0.425E0 / 27.21E0, sqrt(1060.)); };
		Applyfx(HPsi, Psi, grid.GetX(), f);
	}

	double ProductProjector(double x) {
		double value = 1.;

//		double alpha = 5.;
//		value = 1./(1+exp(-alpha*x));

		if ((x < 0)) { value = 0.; }

		return value;
	}

	//! Generell Commutator C=[A,B] for two SPOs acting on the same mode
	/*!
	 \param A The first SPO
	 \param B The second SPO
	 \param Psi the SPFs
	 \result The result is [A,B]Psi.
	 */
	Tensorcd Commutator(const PrimitiveBasis& grid, const Tensorcd& Psi,
		function<Tensorcd(const PrimitiveBasis&, const Tensorcd& Psi)> A,
		function<Tensorcd(const PrimitiveBasis&, const Tensorcd& Psi)> B) {

		// Apply AB
		Tensorcd APsi = A(grid, Psi);
		Tensorcd BAPsi = B(grid, APsi);

		// Apply BA
		Tensorcd BPsi = B(grid, Psi);
		Tensorcd ABPsi = A(grid, BPsi);

		// C = [A, B] = AB-BA
		Tensorcd CPsi(Psi);
		const TensorDim& tdim = CPsi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); n++)
			for (size_t j = 0; j < tdim.getdimpart(); j++)
				CPsi(j, n) = ABPsi(j, n) - BAPsi(j, n);

		return CPsi;
	}

	//! Apply the dividing surface, that is given by the function ProductProjector
	Tensorcd Applyh(const PrimitiveBasis& grid, const Tensorcd& Psi) {
		const Vectord& x = grid.GetX();
		Tensorcd hPsi(Psi);
		Applyfx(hPsi, Psi, x, ProductProjector);
		return hPsi;
	}

	//! Apply the Flux-Operator to a Tensor
	void ApplyFlux(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi) {
		function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;

		// 1. Apply the commutator [T, h]
		HPsi = Commutator(grid, Psi, kin, Applyh);

		// 2. Apply the factor -i
//		complex<double> im(0., 1.);
//		HPsi*=im;
	}
}

//! Initialize the scattering-operator
void scattering::SpecialInitialize(const mctdhBasis& basis) {
	using namespace scatteringoperators;

	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::ApplyX2;

	// part 1: kinetic energy
	{
		MultiParticleOperator M;
		// p�/2
		M.push_back(kin, 0);
		push_back(M, 1.);
	}

	// part 2: eckhard barrier
	{
		MultiParticleOperator M;
		M.push_back(V1, 0);
		push_back(M, 1.);
	}

/*	// Second dimension
	{
		// p�/2
		MultiParticleOperator M;
		M.push_back(kin, 1);
		push_back(M, 1.);
	}

	// 
	{
		MultiParticleOperator M;
		M.push_back(x, 1);
		push_back(M, 0.5E0 * pow(0.2E0 / 27.21E0, 2));
	}
 */
}

