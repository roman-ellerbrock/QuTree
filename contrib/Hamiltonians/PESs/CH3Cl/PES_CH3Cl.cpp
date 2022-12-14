#include "PES_CH3Cl.h"

extern "C"
{
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	void __stdcall POTENTIALCH3CL(double* x, double* w, int* part);
#endif
#ifdef __linux__
	void potentialch3cl_(double* x, double* w, int* part);
#endif
}

template<size_t i2, size_t j2>
complex<double> project(size_t i1, size_t j1)
{
	if ((i1 == i2) && (j1 == j2)){
		return 1.;
	} else if ((i1 == j2) && (j1 == i2)) {
		return 1.;
	} else {
		return 0.;
	}
}

template<size_t l, size_t m>
void ProjectElectronic(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi) {
	assert(HPsi.Dim() == Psi.Dim());
	HPsi.Zero();
	for (size_t n = 0; n < Psi.Dim().getntensor(); ++n)
		for (size_t i = 0; i < Psi.Dim().getdimpart(); ++i)
			for (size_t j = 0; j < Psi.Dim().getdimpart(); ++j)
				HPsi(i, n) += project<l, m>(i, j) * Psi(j, n);
}

template <size_t j, size_t i>
MultiParticleOperator DiabaticOperator(size_t n_el_states, size_t k_elec)
{
	MultiParticleOperator M;
	M.push_back(ProjectElectronic<j, i>, k_elec);
	size_t part = j * n_el_states + i;
	PotentialOperator V(k_elec, part);
	M.SetV(V);

	return M;
}

vector<MultiParticleOperator> PES_CH3Cl::Build(const mctdhBasis& basis)
{
	// mode of the electronic basis
	size_t k_elec = basis.nLeaves() - 1;
	vector<MultiParticleOperator> diab_op;
	size_t nstates = 3;

	// Diagonals
	diab_op.push_back(DiabaticOperator<0, 0>(nstates, k_elec));
	diab_op.push_back(DiabaticOperator<1, 1>(nstates, k_elec));
	diab_op.push_back(DiabaticOperator<2, 2>(nstates, k_elec));

	// Off-Diagonals
	diab_op.push_back(DiabaticOperator<0, 1>(nstates, k_elec));
	diab_op.push_back(DiabaticOperator<0, 2>(nstates, k_elec));
	diab_op.push_back(DiabaticOperator<1, 2>(nstates, k_elec));

	return diab_op;
}

double PES_CH3Cl::Evaluate(const Vectord& Xv, size_t part)
{
	double w = 0;
	// Copying and Converting to fit to fortran format
	Vectord Xv2(Xv);
	assert(Xv.Dim()==15);
	Xv2(0) = Xv(12);
	Xv2(1) = Xv(13);
	Xv2(2) = Xv(14);
	Xv2(12) = Xv(0);
	Xv2(13) = Xv(1);
	Xv2(14) = Xv(2);
	for (size_t i = 0; i < Xv2.Dim(); i++)
		Xv2(i) *= 0.5291777;
	int part_ = (int) part;
//	assert(part_ == 0);
#ifdef _MSC_VER
	POTENTIALCH3CL(&Xv2(0), &w, &part_);
#endif
#ifdef __linux__
	potentialch3cl_(&Xv2(0), &w, &part_);
#endif
//	for (size_t i = 0; i < Xv2.Dim(); i++)
//		cout << Xv2(i) << "\t";
//	cout << w << endl;
	return w;
}

