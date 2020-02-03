#pragma once
#include "Matrix.h"
#include "Tensor.h"

template <typename T>
class DoubleParticleOperator:
	public Matrix<T>
{
public:
	DoubleParticleOperator(size_t nl, size_t nr, size_t k1, size_t k2)
		:Matrix<T>(nl*nl,nr*nr), mode1_(k1), mode12_(k2), nl_(nl), nr_(nr)
	{
		assert(nl_ > 0);
		assert(nr_ > 0);
	}

	T operator()(size_t il, size_t jl, size_t ir, size_t jr)const
	{
		size_t L = nl_ * jl + il;
		size_t R = nr_ * jr + ir;
		return Matrix<T>::operator()(L, R);
	}

	T& operator()(size_t il, size_t jl, size_t ir, size_t jr)
	{
		size_t L = nl_ * jl + il;
		size_t R = nr_ * jr + ir;
		assert(L < ProdDimL());
		assert(R < ProdDimR());
		return Matrix<T>::operator()(L, R);
	}

	size_t DimL()const { return nl_; }
	size_t DimR()const { return nr_; }

	size_t ProdDimL()const { return nl_*nl_; }
	size_t ProdDimR()const { return nr_*nr_; }

	size_t Mode1()const { return mode1_; }
	size_t Mode2()const { return mode12_; }

protected:
	size_t nl_, nr_;
	size_t mode1_, mode12_;

};

typedef DoubleParticleOperator<complex<double>> DPOcd;
typedef DoubleParticleOperator<double> DPOd;

template <typename T>
DoubleParticleOperator<T> TensorDoubleHoleSelect(const Tensor<T>& A,
	const Tensor<T>& B, size_t k1, size_t k2, size_t BathIdx)
{
	
	const TensorDim& tdim = A.Dim();
	const TensorDim& tdim_b = B.Dim();
	assert(tdim == tdim_b);
	size_t active1 = tdim.Active(k1);
	size_t active2 = tdim.Active(k2);

	size_t numgrid = tdim.getdimpart();
	size_t before = tdim.Before(k1);
	size_t behind = tdim.Behind(k2);
	size_t middle = numgrid / (before * behind * active1 * active2);

	// Split DoubleBathIdx into bef, mid, beh
	// @TODO: Should be ok, but double-check this
	size_t before2 = tdim.Before(k2);
	size_t beh = BathIdx / before2;
	BathIdx = BathIdx % before2;
	size_t mid = BathIdx / middle;
	size_t bef = BathIdx % middle;
	cout << "bef, mid, beh = " << bef << " " << mid << " " << beh << endl;

	DoubleParticleOperator<T> zeta(active1, active2, k1, k2);

	size_t ntensor = tdim.getntensor();
						for (size_t m = 0; m < ntensor; m++)
						{
					for (size_t i2 = 0; i2 < active2; i2++)
					{
				for (size_t j2 = 0; j2 < active2; j2++)
				{
			for (size_t i = 0; i < active1; i++)
			{
		for (size_t j = 0; j < active1; j++)
		{
	zeta(i, j, i2, j2)+= 0;
//		conj(A(bef,i,mid,j,beh,k1,k2,m)) * B(bef,i2,mid,j2,beh,k1,k2,m);
		}
			}
				}
					}
						}

	return zeta;
}

