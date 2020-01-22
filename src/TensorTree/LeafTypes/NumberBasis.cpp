#include "NumberBasis.h"



NumberBasis::NumberBasis(int dim_,bool fermion_)
	:PrimitiveBasis()
{
	dim=dim_;
	fermion=fermion_;

	//check basis size for fermions
	if(fermion) assert(dim <= 2);
}

NumberBasis::~NumberBasis()
{
}

void NumberBasis::Initialize(double occ,
                             double oSQR,
                             double minimalOccupation,
                             double dummy3)
{  
	// save all parameters
	startocc = (int)occ;
	osqrtree = (int)oSQR;
	minOcc = (int)minimalOccupation;
	
	if(minOcc > startocc)
	{
		cout << "Error: Minimal occupation must not be bigger than\n"
		     << "the initial occupation.\n";
		assert(0);
	}
}

void NumberBasis::InitSPF(Tensorcd& phi)const
{
	TensorDim tdim(phi.Dim());
	int nstates = tdim.getntensor();

	if(fermion && nstates > 2)
	{
		cout << "Error: there can be only two fermionic states\n"
		     << "per mode.\n";
		assert(0);
	}
	if(fermion && minOcc > 2)
	{
		cout << "Error: For fermions the minimal Occupation must be\n"
		     << "Zero.\n";
		assert(0);
	}

	// soft check for bottom layer_
	assert(tdim.F() == 1);
	assert(tdim.getdimpart() == dim);
	assert(startocc - minOcc < dim);
	
	// set ground state wf
	for (int i = 0; i < dim; i++)
	{
		phi(i, 0) = 1.e-7;
	}
	phi(startocc - minOcc , 0) = 1.0;
  
	// excitations
	int count = 0;
	for (int n = 1; n < nstates; n++)
	{
		for (int i = 0; i < dim; i++)
		{
			phi(i,n) = 0.0;
		}
		if(count == startocc - minOcc) count++;
		phi(count,n) = 1.0;
		count++;
	}
	// orthonormalize
	GramSchmidt(phi);
}

Tensorcd NumberBasis::ToGrid(const Tensorcd& phi)const
{
	// soft check that its really a bottom-layer_ tensor
	TensorDim tdim = phi.Dim();
	assert(tdim.F() == 1);

	return phi;
}

Tensorcd NumberBasis::FromGrid(const Tensorcd& phi)const
{
	// soft check that its really a bottom-layer_ tensor
	TensorDim tdim = phi.Dim();
	assert(tdim.F() == 1);

	return phi;
}

Tensorcd NumberBasis::ApplyKin(const Tensorcd& phi)const
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer_ tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	int nstates = tdim.getntensor();
	int active = tdim.getdimpart();
	
	assert(active == dim);

	for (int n = 0; n < nstates; n++)
	{
		for (int i = 0; i < active; i++)
		{
			psi(i, n) = (1.0*(i + minOcc))*phi(i, n);
		}
	}
	return psi;
}

Tensorcd NumberBasis::ApplyP(const Tensorcd& phi)const
{
	const TensorDim& tdim = phi.Dim();
	Tensorcd psi(tdim, false);

	size_t prim = tdim.getdimpart();
	size_t states =  tdim.getntensor();

	for (size_t n = 0; n < states; n++)
	{
		psi[n*prim]=0.0;
		for (size_t i = 1; i < prim; i++)
		{
			psi[n*prim + i] = sqrt(1.*(i + minOcc))*phi[n*prim + i - 1];
		}
	}
	return psi;
}


// Apply primitive x for several single particle functions
Tensorcd NumberBasis::applyX(const Tensorcd& phi)const
{
	const TensorDim& tdim = phi.Dim();
	Tensorcd psi(tdim, false);

	size_t prim = tdim.getdimpart();
	size_t states =  tdim.getntensor();
	
	for (int n = 0; n < states; n++)
	{
		psi[n * prim + prim - 1] = 0.0;
		for (int i = 1; i < prim; i++)
		{
			psi[n*prim + i - 1] = sqrt(1.*(i + minOcc))*phi[n*prim + i];
		}
	}
	return psi;
}
