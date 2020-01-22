#include "DVRBasis.h"

DVRBasis::DVRBasis(int dim_)
  :dim(dim_), trafo(dim_, 0), x(dim_), kin(dim_, 0), p(dim_, 0)
  {
  }

  Tensorcd DVRBasis::applyX(const Tensorcd& phi)const
  {
  	TensorDim tdim = phi.Dim();
  	// check that its really a bottom-layer_ tensor
  	assert(tdim.F() == 1);

  	Tensorcd psi(phi.Dim());

  	int active = tdim.getdimpart();
  	assert(active == dim);

  //	psi = kin*phi;
  // @TODO: rewrite this code as a matrix*Tensor routine
//  	#pragma omp for
  	for (int n = 0; n < tdim.getntensor(); n++)
  	{
  		for (int i = 0; i < active; i++)
  		{
  				psi(i, n) = x(i)*phi(i, n);
  		}
  	}
  	return psi;
  }

  Tensorcd DVRBasis::ApplyX2(const Tensorcd& phi)const
  {
    TensorDim tdim = phi.Dim();
    // check that its really a bottom-layer_ tensor
    assert(tdim.F() == 1);

    Tensorcd psi(phi.Dim());

    int active = tdim.getdimpart();
    assert(active == dim);

  //	psi = kin*phi;
  // @TODO: rewrite this code as a matrix*Tensor routine
//    #pragma omp for
    for (int n = 0; n < tdim.getntensor(); n++)
    {
      for (int i = 0; i < active; i++)
      {
          psi(i, n) = x(i)*x(i)*phi(i, n);
      }
    }
    return psi;
  }

  Tensorcd DVRBasis::ApplyP(const Tensorcd& phi)const
  {
    // soft check that its really a bottom-layer_ tensor
    TensorDim tdim = phi.Dim();
    assert(tdim.F() == 1);

    return p*phi;
  }

  Tensorcd DVRBasis::ApplyKin(const Tensorcd& phi)const
  {
  	// soft check that its really a bottom-layer_ tensor
  	TensorDim tdim = phi.Dim();
  	assert(tdim.F() == 1);

  	return kin*phi;
  }


  Tensorcd DVRBasis::ToGrid(const Tensorcd& phi)const
  {
  	// soft check that its really a bottom-layer_ tensor
  	TensorDim tdim = phi.Dim();
  	assert(tdim.F() == 1);

  	return multATB(trafo, phi);
  }

  Tensorcd DVRBasis::FromGrid(const Tensorcd& phi)const
  {
  	// soft check that its really a bottom-layer_ tensor
  	TensorDim tdim = phi.Dim();
  	assert(tdim.F() == 1);

  	return trafo*phi;
  }
