#include "LinearVibronicCoupling.h"

LinearVibronicCoupling::LinearVibronicCoupling()
{
}


LinearVibronicCoupling::~LinearVibronicCoupling()
{
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
//	for (size_t n = 0; n < Psi.Dim().getntensor(); ++n)
//		for (size_t i = 0; i < Psi.Dim().getdimpart(); ++i)
//			HPsi(i, n) += Psi(i, n);

	for (size_t n = 0; n < Psi.Dim().getntensor(); ++n)
		for (size_t i = 0; i < Psi.Dim().getdimpart(); ++i)
			for (size_t j = 0; j < Psi.Dim().getdimpart(); ++j)
			{
//				cout << "i, j= " << i << " " << j << " " << project<l,m>(i,j)<< endl;
				HPsi(i, n) += project<l, m>(i, j) * Psi(j, n);
			}
//	Psi.print();
//	HPsi.print();
//	getchar();
}

template <size_t j, size_t i>
MultiParticleOperator DiabaticOperator(size_t n_el_states, size_t k_elec)
{
	MultiParticleOperator M;
	M.push_back(ProjectElectronic<j, i>, k_elec);
	size_t part = j * n_el_states + i;
	PotentialOperator V(k_elec, part);
	M.SetV(V);

/*
	cout << "i = " << i << " j = " << j << endl;
	cout << "part = " << part << endl;
	cout << "k_elec = " << k_elec << endl;
	cout << endl;
	getchar();
	*/
	return M;
}

vector<MultiParticleOperator> LinearVibronicCoupling::Build(const mctdhBasis& basis)
{
	// mode of the electronic basis
	size_t k_elec = basis.nLeaves() - 1;
	size_t nstates = 2;
	vector<MultiParticleOperator> diab_op;

	// Diagonals
	diab_op.push_back(DiabaticOperator<0, 0>(nstates, k_elec));
	diab_op.push_back(DiabaticOperator<1, 1>(nstates, k_elec));

	// Off-Diagonal
	diab_op.push_back(DiabaticOperator<0, 1>(nstates, k_elec));

	return diab_op;
}

double LinearVibronicCoupling::Evaluate(const Vectord & Xv, size_t part)
{
	double a = 1.;
	double v = 0;

	if (part == 0)
	{
		v += 0.5*Xv(0)*Xv(0);
		v += 0.5*Xv(1)*Xv(1);
		v += a*Xv(0);
	}
	else if (part == 3)
	{
		v += 0.5*Xv(0)*Xv(0);
		v += 0.5*Xv(1)*Xv(1);
		v -= a*Xv(0);
	}
	else if (part == 1)
	{
		v += a*Xv(1);
	}
	else
	{
		cerr << "Wrong part number." << endl;
		exit(1);
	}

	return v;
}

