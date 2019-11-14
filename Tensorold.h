#pragma once
#include "TensorDim.h"
#include "Matrix.h"

class mctdhBlockWavefunction;

template <typename T>
class Tensor
{
public:
	//////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////
	friend mctdhBlockWavefunction;

	Tensor()
		:coeffs(new T[1])
	{}

	Tensor(const TensorDim dim_, const bool InitZero=true)
		:dim(dim_), coeffs(new T[dim_.getdimtot()])
	{
		if (InitZero) { Zero(); }
	}

	// Copy constructor
	Tensor(const Tensor& old)
		:Tensor(old.dim, false)
	{
		#pragma omp for
		for (int i = 0; i<dim.getdimtot(); i++)
		{
			coeffs[i] = old.coeffs[i];
		}
	}

	// Copy-Multyply constructor
	Tensor(const Tensor& old, T factor)
		:Tensor(old.dim, false)
	{
		#pragma omp for
		for (int i = 0; i<dim.getdimtot(); i++)
		{
			coeffs[i] = old.coeffs[i]*factor;
		}
	}

	// Move constructor
	Tensor(Tensor&& old)
		:dim(old.dim), coeffs(old.coeffs)
	{
		old.coeffs = nullptr;
	}

	// Copy Assignment Operator
	Tensor& operator=(const Tensor& old)
	{
		Tensor tmp(old);
		*this = move(tmp);
		return *this;
	}
	
	// Move Assignment Operator
	Tensor& operator=(Tensor&& old)
	{
		delete[] coeffs;
		dim = old.dim;
		coeffs = old.coeffs;
		old.coeffs = nullptr;
		return *this;
	}


	// Destructor
	~Tensor()
	{
		delete[] coeffs;
	}

	//////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////
	inline T& operator()(const int i)const
	{
		int dimtot = dim.getdimtot();
		assert(i >= 0 && i < dimtot);
		return coeffs[i];
	}

	inline T& operator()(const int i)
	{
		int dimtot = dim.getdimtot();
		assert(i >= 0 && i < dimtot);
		return coeffs[i];
	}

	inline T& operator()(const int i,const int n)const
	{
		int dimpart = dim.getdimpart();
		assert(i >= 0 && i < dimpart);
		assert(n >= 0 && n < dim.getntensor());
		return coeffs[n*dimpart + i];
	}

	inline T& operator()(const int i,const int n)
	{
		int dimpart = dim.getdimpart();
		assert(i >= 0 && i < dimpart);
		assert(n >= 0 && n < dim.getntensor());
		return coeffs[n*dimpart + i];
	}

	inline T& operator()(const int i, const int j, const int k, const int f, const int n)
	{
		assert(i >= 0);
		assert(j >= 0);
		assert(k >= 0);
		assert(f >= 0);
		assert(n >= 0);
		int a = dim.Before(f);
		int b = dim.Active(f);
		int dimpart = dim.getdimpart();
		int idx = n*dimpart + k*a*b + j*a + i;
		assert(i < a);
		assert(j < b);
		assert(n < dim.getntensor());
		assert(f < dim.F());
		assert(k < dim.Behind(f));
		return coeffs[idx];
	}

	inline T& operator()(const int i, const int j, const int k, const int f, const int n)const
	{
		assert(i >= 0);
		assert(j >= 0);
		assert(k >= 0);
		assert(f >= 0);
		assert(n >= 0);
		int a = dim.Before(f);
		int b = dim.Active(f);
		int dimpart = dim.getdimpart();
		int idx = n*dimpart + k*a*b + j*a + i;
		assert(i < a);
		assert(j < b);
		assert(n < dim.getntensor());
		assert(f < dim.F());
		assert(k < dim.Behind(f));
		return coeffs[idx];
	}

	inline T& operator[](const int idx)const
	{
		// Fast bracket operator
		return coeffs[idx];
	}

	inline T& operator[](const int idx)
	{
		// Fast bracket operator
		return coeffs[idx];
	}

	void WriteRaw(ostream& os)const
	{
		const char record = 255;
		for (int i = 0; i < dim.getdimtot(); i++)
		{
			os << operator()(i);
			for (int j=0; j<8; j++)
				os << record;
		}
	}

	void ReadRaw(istream& is)
	{
		char record;
		for (int i = 0; i < dim.getdimtot(); i++)
		{
			is >> operator()(i);
			for (int j=0; j<8; j++)
				is >> record;
		}
	}

	void Write(ostream& os)const
	{
		for (int i = 0; i < dim.getdimtot(); i++)
			os << operator()(i);
	}

	void WriteBin(ofstream& os)const
	{
		// Verification
		os.write("TENS", 4);

		// Write the TensorDim
		dim.WriteBin(os);

		// Write the size
		int32_t size = sizeof(T);
		os.write((char*)&size, sizeof(size));

		// Write the Coefficients
		for (int i = 0; i < dim.getdimtot(); i++)
		{
			T Coeff_now = operator()(i);
			os.write((char*) &Coeff_now, size);
		}
	}

	void ReadBin(ifstream& is)
	{
		// Check if binary string contains a Tensor
		char check[5];
		is.read(check, 4);
		string s_check(check, 4);
		string s_key("TENS");
		assert(s_key == s_check);

		// Read the TensorDim
		TensorDim newtdim;
		newtdim.ReadDim(is);

		// Resize the Tensor
		(*this) = Tensor<T>(newtdim, false);

		// Read the size
		int32_t size;
		is.read((char*)&size, sizeof(size));
		assert(size == sizeof(T));

		// Read the coefficients

		for (int i = 0; i < dim.getdimtot(); i++)
		{
			T Coeff_now;
			is.read((char*)&Coeff_now, size);
			operator()(i) = Coeff_now;
		}
	}

	void Read(istream& os)
	{
		for (int i = 0; i < dim.getdimtot(); i++)
			os >> operator()(i);
	}

	friend Tensor operator+(const Tensor& A, const Tensor& B)
	{
		assert(A.Dim().getdimtot() == B.Dim().getdimtot());
		Tensor C(A.Dim());
		#pragma omp for
		for (int i = 0; i < A.Dim().getdimtot(); i++)
		{
			C(i) = A(i) + B(i);
		}
		return C;
	}

	void operator+=(const Tensor& A)
	{
		assert(A.Dim().getdimtot() == Dim().getdimtot());
		#pragma omp for
		for (int i = 0; i < A.Dim().getdimtot(); i++)
		{
			(*this)(i) += A(i);
		}
	}

	friend Tensor coeffprod(const Tensor& A, const Tensor& B)
	{
		assert(A.Dim().getdimtot() == B.Dim().getdimtot());
		Tensor C(A.Dim());
		for (int i = 0; i < A.Dim().getdimtot(); i++)
		{
			C(i) = A(i)*B(i);
		}
		return C;
	}

	void operator*=(T a)
	{
		#pragma omp for
		for (int i = 0; i < Dim().getdimtot(); i++)
		{
			operator()(i) = a*operator()(i);
		}
	}

	friend Tensor operator-(const Tensor& A, const Tensor& B)
	{
		Tensor C(A.Dim());
		#pragma omp for
		for (int i = 0; i < A.Dim().getdimtot(); i++)
		{
			C(i) = A(i) - B(i);
		}
		return C;
	}

	Tensor<T> AdjustDimensions(const TensorDim& newTDim)
	{
		// Increase the dimensions of the Tensor from old TensorDim
		// to new TensorDim 

		assert(newTDim.F() == dim.F());

		// Increase the active modes
		Tensor<T> Acoeff(*this);
		for (int k = 0; k < dim.F(); k++)
		{
			int act = newTDim.Active(k);
			Acoeff = Acoeff.AdjustActiveDim(act, k);
		}

		// Increase the number of Tensors
		int ntens = newTDim.getntensor();
		Acoeff = Acoeff.AdjustStateDim(ntens);

		return Acoeff;

	}

	Tensor<T> AdjustActiveDim(int active, int mode)
	{
		// Adjust the active dimension in the coordinate "mode".
		// If the new active is smaller, the norm of the tensors is
		// not conserved.

		assert(0 <= mode);
		assert(mode < dim.F());
		assert(active > 0);

		// Create a new Tensor with the adjusted dim
		vector<int> dimlist = dim.getdimlist();
		dimlist[mode] = active;
		int ntensor = dim.getntensor();
		TensorDim newTDim(dimlist, ntensor);
		Tensor<T> newT(newTDim);

		// Copy the coefficients
		int before = dim.Before(mode);
		int behind = dim.Behind(mode);
		int minactive = min(active, dim.Active(mode));
		for (int n = 0; n < ntensor; n++)
		{
			for (int l = 0; l < behind; l++)
			{
				for (int j = 0; j < minactive; j++)
				{
					for (int i = 0; i < before; i++)
					{
						newT(i, j, l, mode, n) = operator()(i, j, l, mode, n);
					}
				}
			}
		}
		return newT;
	}

	// Adjust the size of Tensor 
	Tensor<T> AdjustStateDim(int n)
	{
		// Returns a new tensor with n (>=ntensor) tensors
		// The new tensors are all set to zero

		// Create a new TensorDim with the new size
		vector<int> dimlist = dim.getdimlist();
		TensorDim newTDim(dimlist, n);
		Tensor<T> newTensor(newTDim);

		// Copy the coefficients
		int ntensor = dim.getntensor();
		int dimpart = dim.getdimpart();
		int ntensmax = min(n, ntensor);
		for (int m = 0; m < ntensmax; m++)
		{
			for (int i = 0; i < dimpart; i++)
			{
				newTensor(i, m) = operator()(i, m);
			}
		}
		return newTensor;
	}

	// Getter for Dim
	const TensorDim& Dim()const { return dim; }
	TensorDim& Dim() { return dim; }

	void Zero()
	{
		#pragma omp for
		for (int i = 0; i < dim.getdimtot(); i++)
			coeffs[i] = 0;
	}

	void print(ostream& os=cout)const
	{
		for (int n = 0; n < dim.getntensor(); n++)
		{
			for (int i = 0; i < dim.getdimpart(); i++)
				os << (*this)(i, n) << " ";
			os << endl;
		}
		os << endl;
	}

	T singleDotProduct(const Tensor& A, int n, int m)const
	{
		T result = 0;
		for (int i = 0; i < A.Dim().getdimpart(); i++)
		{
			result += conj(operator()(i, n))*A(i, m);
		}
		return result;
	}

	Matrixcd DotProduct(const Tensor<T>& A)const
	{
		TensorDim tdima(A.Dim());
		// Every tensor can have different amount of states but same dimpart
		assert(tdima == dim);
	
		int nmax = tdima.getntensor();
		int mmax = dim.getntensor();
		int npart = dim.getdimpart();
	
		Matrixcd S(nmax, mmax);
		#pragma omp for
		for (int n = 0; n < nmax; n++)
		{
			for (int m = 0; m < mmax; m++)
			{
				for (int i = 0; i < npart; i++)
				{
					S(m, n) += conj(operator()(i, m))*A(i, n);
				}
			}
		}
		return S;
	}

protected:
	TensorDim dim;
	T* coeffs;
};


template<typename T>
T SingleDotProd(const Tensor<T>& A, const Tensor<T>& B, int n, int m)
{
	TensorDim tdima(A.Dim());
	TensorDim tdimb(B.Dim());
	// Every tensor can have different amount of states but same dimpart
	assert(tdima.getdimpart() == tdimb.getdimpart());

	int nmax = tdima.getntensor();
	int mmax = tdimb.getntensor();
	int npart = tdima.getdimpart();
	assert(n >= 0 && n < nmax);
	assert(m >= 0 && m < mmax);

	T result = 0;
	#pragma omp for
	for (int i = 0; i < npart; i++)
	{
		result += conj(A(i, n))*B(i, m);
	}
	return result;
}

template<typename T>
void TensorHoleProduct(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	int before, int active, int behind)
{
	// Variables for precalculation of indices
	int actbef = active*before;
	int Sidx = 0;
	int Aidx = 0;
	int Bidx = 0;
	int ipreidx = 0;
	int jpreidx = 0;
	int npreidx = 0;

	#pragma omp for
	for (int n = 0; n < behind; n++)
	{
		npreidx = n*actbef;
		for (int j = 0; j < active; j++)
		{
			jpreidx = npreidx + j*before;
			for (int i = 0; i < active; i++)
			{
				// S(i, j)
				Sidx = j*active + i;
				ipreidx = npreidx + i*before;
				for (int l = 0; l < before; l++)
				{
					// A(l, i, n)
					Aidx = ipreidx + l;
					// B(l, j, n)
					Bidx = jpreidx + l;
					S[Sidx] += conj(A[Aidx]) * B[Bidx];
				}
			}
		}
	}

}

template<typename T>
Matrix<T> StateAveragedHoleProduct(const Tensor<T>& A, const Tensor<T>& B, int k)
{
	TensorDim tdim(A.Dim());
	// check wether tensordims are equal
	assert(tdim == B.Dim());

	int nstates = tdim.getntensor();
	int active = tdim.Active(k);
	int before = tdim.Before(k);
	int behind = tdim.Behind(k)*nstates;
	Matrix<T> S(active, active);

	TensorHoleProduct(S, A, B, before, active, behind);

	return S;
}

template<typename T>
Matrix<T> OldStateAveragedHoleProduct(const Tensor<T>& A, const Tensor<T>& B, int k)
{
	TensorDim tdim(A.Dim());
	// check wether tensordims are equal
	assert(tdim == B.Dim());

	int nstates = tdim.getntensor();
	int active = tdim.Active(k);
	int before = tdim.Before(k);
	int behind = tdim.Behind(k);
	Matrix<T> S(active, active);

	for (int n = 0; n < nstates; n++)
	{
		for (int beh = 0; beh < behind; beh++)
		{
			for (int i = 0; i < active; i++)
			{
				for (int j = 0; j < active; j++)
				{
					for (int bef = 0; bef < before; bef++)
					{
						S(i, j) += conj(A(bef, i, beh, k, n))*B(bef, j, beh, k, n);
					}
				}
			}
		}
	}

	// divide by number of states
//	for (int j = 0; j < active; j++)
//		for (int i = 0; i < active; i++)
//			S(i, j) /= nstates;

	return S;
}


template<typename T>
void GramSchmidt(Tensor<T>& A)
{
	// @TODO: Fill in auto-refill

	// control parameters
	int maxiter = 15;
	double conver = 1e-12;
	double errorconver = 1e-9;

	TensorDim tdim(A.Dim());
	int ntensor = tdim.getntensor();
	int dimpart = tdim.getdimpart();

	for (int n = 0; n < ntensor; n++)
	{
		int iter = 0;
		double accumoverlap = 1.;
		// orthogonalize on all previous ones and then normalize
		while ((accumoverlap > conver) && (iter < maxiter))
		{
			iter++;
			accumoverlap = 0;
			for (int m = 0; m < n; m++)
			{
				// orthogonalize
				T overlap = SingleDotProd(A, A, m, n);
				accumoverlap += abs(overlap);
				for (int i = 0; i < dimpart; i++)
				{
					A(i, n) -= overlap*A(i, m);
				}
			}

			// renormalize
			T norm = SingleDotProd(A, A, n, n);
			if (abs(norm) != 0)
			{
				norm = sqrt(real(norm));
				for (int i = 0; i < dimpart; i++)
				{
					A(i, n) /= norm;
				}
			}
		}
		// Error message
		if (accumoverlap >= errorconver) {
			cout << "Error: No orthogonality in Gram-Schmidt" << endl;
			cout << "Error measurement: " << conver << endl;
			assert(0);
		}
	}
}

template <typename T, typename U>
void mattensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	int before, int active, int behind, bool zero=true)
{
	// Null the result tensor if flag is set to "true"
	if (zero) { C.Zero(); }

	// Variables to Precompute index values
	int actbef = active*before;
	int Aidx = 0;
	int Bidx = 0;
	int Cidx = 0;
	int kpreidx = 0;
	int lpreidx = 0;
	int jpreidx = 0;
	int lactive = 0;

	if (before == 1)
	{
		#pragma omp for
		for (int k = 0; k < behind; k++)
		{
			kpreidx = k*actbef;
			for (int l = 0; l < active; l++)
			{
				Bidx = l + kpreidx;
				for (int j = 0; j < active; j++)
				{
					Cidx = j + kpreidx;
					Aidx = l*active + j;
					C[Cidx] += A[Aidx] * B[Bidx];
				}
			}
		}
	}
	else
	{
		#pragma omp for
		for (int k = 0; k < behind; k++)
		{
			kpreidx = k*actbef;
			for (int l = 0; l < active; l++)
			{
				lpreidx = l*before + kpreidx;
				lactive = l*active;
				for (int j = 0; j < active; j++)
				{
					Aidx = lactive + j;
					jpreidx = j*before + kpreidx;
					for (int i = 0; i < before; i++)
					{
						Cidx = jpreidx + i;
						Bidx = lpreidx + i;
						C[Cidx] += A[Aidx] * B[Bidx];
					}
				}
			}
		}
	}
}

template <typename T, typename U>
void Tmattensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	int before, int active, int behind, bool zero = true)
{
	// Null the result tensor if flag is set to "true"
	if (zero) { C.Zero(); }

	/*	int actbef = active*before;
		int Cidx = 0;
		int Bidx = 0;
		int Aidx = 0;
		int kpreidx = 0;
		#pragma omp for
		for (int k = 0; k < behind; k++)
		{
			kpreidx = k*actbef;
			for (int l = 0; l < active; l++)
				for (int j = 0; j < active; j++)
					for (int i = 0; i < before; i++)
					{
						Cidx = kpreidx + j*before + i;
						Bidx = kpreidx + l*before + i;
						Aidx = j*active + l;
						C[Cidx] += conj(A[Aidx]) * B[Bidx];
					}
		}
		*/
	int actbef = active*before;
	int Cidx = 0;
	int Bidx = 0;
	int Aidx = 0;
	int kpreidx = 0;
	int lpreidx = 0;
	int jpreidx = 0;

	if (before == 1)
	{
	#pragma omp for
		for (int k = 0; k < behind; k++)
		{
			kpreidx = k*actbef;
			for (int l = 0; l < active; l++)
			{
				Bidx = l + kpreidx;
				for (int j = 0; j < active; j++)
				{
					Cidx = j + kpreidx;
					Aidx = j*active + l;
					C[Cidx] += conj(A[Aidx]) * B[Bidx];
				}
			}
		}
	}
	else
	{
		#pragma omp for
		for (int k = 0; k < behind; k++)
		{
			kpreidx = k*actbef;
			for (int l = 0; l < active; l++)
			{
				lpreidx = l*before + kpreidx;
				for (int j = 0; j < active; j++)
				{
					Aidx = j*active + l;
					jpreidx = j*before + kpreidx;
					for (int i = 0; i < before; i++)
					{
						Cidx = jpreidx + i;
						Bidx = lpreidx + i;
						C[Cidx] += conj(A[Aidx]) * B[Bidx];
					}
				}
			}
		}
	}
}

template <typename T, typename U>
void multAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, int mode, bool zero=true)
{
	TensorDim tdim(B.Dim());
	assert(mode < tdim.F());
	assert(mode >= 0);
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim1() == B.Dim().Active(mode));

	int behind = tdim.Behind(mode)*tdim.getntensor();
	int active = tdim.Active(mode);
	int before = tdim.Before(mode);
	mattensor(C, A, B, before, active, behind, zero);
}

template <typename T, typename U>
Tensor<T> multAB(const Matrix<U>& A, const Tensor<T>& B, int mode)
{
	TensorDim tdim(B.Dim());
	assert(mode < tdim.F());
	assert(mode >= 0);
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim1() == B.Dim().Active(mode));

	Tensor<T> C(tdim);
	int behind = tdim.Behind(mode)*tdim.getntensor();
	int active = tdim.Active(mode);
	int before = tdim.Before(mode);
	mattensor(C, A, B, before, active, behind, false);
	return C;
}

template <typename T, typename U>
Tensor<T> OldmultAB(const Matrix<U>& A, const Tensor<T>& B, int mode)
{
	TensorDim tdim(B.Dim());
	assert(mode < tdim.F());
	assert(mode >= 0);
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim1() == B.Dim().Active(mode));

	Tensor<T> C(tdim);
	#pragma omp for
	for (int n=0; n<tdim.getntensor(); n++)
		for (int k = 0; k < tdim.Behind(mode); k++)
			for (int l=0; l<tdim.Active(mode); l++)
				for (int j=0; j<tdim.Active(mode); j++)
				for (int i=0; i<tdim.Before(mode); i++)
				{
					C(i, j, k, mode, n) += A(j, l)*B(i, l, k, mode, n);
				}
	return C;
}

template <typename T, typename U>
Tensor<T> multATB(const Matrix<U>& A, const Tensor<T>& B, int mode)
{
	TensorDim tdim(B.Dim());
	assert(mode < tdim.F());
	assert(mode >= 0);
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim1() == B.Dim().Active(mode));

	Tensor<T> C(tdim);
	int behind = tdim.Behind(mode)*tdim.getntensor();
	int active = tdim.Active(mode);
	int before = tdim.Before(mode);
	Tmattensor(C, A, B, before, active, behind, false);
	return C;
	
}

template <typename T, typename U>
void multStateAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, bool zero=true)
{
	const TensorDim& tdim(B.Dim());
	int ntensor = tdim.getntensor();
	int dimpart = tdim.getdimpart();
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim2() == ntensor);

	const int before=dimpart;
	const int active=ntensor;
	const int behind=1;

	mattensor(C, A, B, before, active, behind, zero);
}

template <typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B)
{
	const TensorDim& tdim(B.Dim());
	int ntensor = tdim.getntensor();
	int dimpart = tdim.getdimpart();
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim2() == ntensor);

	Tensor<T> C(tdim);
	int Cidx = 0;
	int Aidx = 0;
	int Bidx = 0;
	int mpreidx = 0;
	int npreidx = 0;
	int apreidx = 0;
	for (int n = 0; n < ntensor; n++)
	{
		npreidx = n*dimpart;
		apreidx = n*ntensor;
		for (int m = 0; m < ntensor; m++)
		{
			mpreidx = m*dimpart;
			Aidx = apreidx + m;
			#pragma omp for
			for (int i = 0; i < dimpart; i++)
			{
				Cidx = mpreidx + i;
				Bidx = npreidx + i;
				C[Cidx] += A[Aidx]*B[Bidx];
			}
		}
	}
	return C;
}

template <typename T, typename U>
Tensor<T> OldmultStateAB(const Matrix<U>& A, const Tensor<T>& B)
{
	TensorDim tdim(B.Dim());
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim2() == B.Dim().getntensor());

	Tensor<T> C(tdim);
	for (int n=0; n<tdim.getntensor(); n++)
		for (int m=0; m<tdim.getntensor(); m++)
			for (int i=0; i<tdim.getdimpart(); i++)
				C(i, m) += A(m, n)*B(i, n);

	return C;
}

template <typename T, typename U>
Tensor<T> multStateArTB(const Matrix<U>& A, const Tensor<T>& B)
{
	TensorDim tdim(B.Dim());
	assert(A.Dim1() == A.Dim2());
	assert(A.Dim2() == B.Dim().getntensor());

	Tensor<T> C(tdim);
	for (int n=0; n<tdim.getntensor(); n++)
		for (int m=0; m<tdim.getntensor(); m++)
			for (int i=0; i<tdim.getdimpart(); i++)
				C(i, m) += A(n, m)*B(i, n);

	return C;
}

//Contracts two tensors A and B over the indeces a and b
//Only the second tensor is allowed to have several staed
template <typename T, typename U>
Tensor<T> TensorProduct(const Matrix<U>& A, const Tensor<T>& B, int a, int b)
{
	TensorDim& dimA = A.Dim();
	TensorDim& dimB = B.Dim();

	assert(dimA.getntensor() == 1);
	assert(a >= 0);
	assert(b >= 0);
	assert(a < dimA.F());
	assert(b < dimB.F());

	vector< int > newdim;
	for(int i = 0; i < dimA.F(); i++)
	{
		newdim.push_back(dimA.Active(i));
	}
	for(int i = 0; i < dimB.F(); i++)
	{
		newdim.push_back(dimB.Active(i));
	}
	
	TensorDim dimC(newdim, dimB.getntensor());
	Tensor<T> C(dimC, false);

	 int beforeA = dimA.Before(a);
	 int beforeB = dimB.Before(b);
	 int afterA = dimA.Behind(a);
	 int afterB = dimB.Behind(b);
	
	//... Multiplication

	return C;
}

template<typename T>
void reverseGramSchmidt(Tensor<T>& A)
{
	// control parameters
	int maxiter = 15;
	double conver = 1e-12;
	double errorconver = 1e-9;

	TensorDim tdim(A.Dim());
	int ntensor = tdim.getntensor();
	int dimpart = tdim.getdimpart();

	for (int n = ntensor-1; n >= 0; n--)
	{
		int iter = 0;
		double accumoverlap = 1.;
		// orthogonalize on all previous ones and then normalize
		while ((accumoverlap > conver) && (iter < maxiter))
		{
			iter++;
			accumoverlap = 0;
			for (int m = ntensor-1; m > n; m--)
			{
				// orthogonalize
				T overlap = SingleDotProd(A, A, m, n);
				accumoverlap += abs(overlap);
				for (int i = 0; i < dimpart; i++)
				{
					A(i, n) -= overlap*A(i, m);
				}
			}

			// renormalize
			T norm = SingleDotProd(A, A, n, n);
			if (abs(norm) != 0)
			{
				norm = sqrt(real(norm));
				for (int i = 0; i < dimpart; i++)
				{
					A(i, n) /= norm;
				}
			}
		}
		// Error message
		if (accumoverlap >= errorconver) {
			cout << "Error: No orthogonality in Gram-Schmidt" << endl;
			cout << "Error measurement: " << conver << endl;
			assert(0);
		}
	}
}

//Projects B on A
template<typename T>
Tensor< complex<double> > Project(const Tensor< complex<double> >& A, 
                                  const Tensor< T >& B)
{
	//calculates the overlap of A with it self
	const Matrix< complex<double> > overlap = A.DotProduct(A);
	
	//invert the overlap
	const Matrix< complex<double> > inverse_operlap = overlap.cInv();
	
	//calculate the scalar product of A and B
	const Matrix< complex<double> > skalarproduct = A.DotProduct(B);

	//multiply the scalar product and the inverse_operlap
	const Matrix< complex<double> > product = inverse_operlap*skalarproduct;

	return multStateArTB(product, A);
}

typedef Tensor<complex<double>> Tensorcd;
typedef Tensor<double> Tensord;
