#pragma once
#include "Core/Vector.h"
#include "Core/stdafx.h"
#include "Core/Tensor.h"
#include <bitset>

template <typename T>
T reverse(T n, size_t b = sizeof(T) * CHAR_BIT)
{
//	assert(b <= std::numeric_limits<T>::digits);

	T rv = 0;

	for (size_t i = 0; i < b; ++i, n >>= 1) {
		rv = (rv << 1) | (n & 0x01);
	}

	return rv;
}



class FFTCooleyTukey
{
public:
	FFTCooleyTukey(){}
	~FFTCooleyTukey(){}

	void Reverse(Vectorcd& A)
	{
		// Perform a bit-reversal on a Vector
		for (int i = 0; i < A.Dim(); i++)
		{
//			unsigned long long b = (unsigned long long) A(i);
			std::bitset<16 * 8> a(100);
		}
	}

	void Init(int dim)
	{
		dim_ = dim;

		logdim_ = 0;
		int rest = dim_;
		while ((rest % 2) == 0)
		{
			logdim_++;
			rest /= 2;
		}
		assert(rest % 2 == 1);

		trafo_ = Vectorcd(dim_);
		trafoback_ = Vectorcd(dim_);
		complex<double> im(0, 1);
		double pi = 3.14159265359;
		for (int i = 0; i < dim_; i++)
			trafo_(i) = exp(im*2.*pi*(1.*i) / (1.*dim_));

		for (int i = 0; i < dim_; i++)
			trafoback_(i) = exp(-im*2.*pi*(1.*i) / (1.*dim_));
	}

	void fft(Vectorcd& A)
	{
		/*
		Vectorcd B(A.Dim());
		int divisor = 2;
		int offset = 1;
		int rest = 0;
		for (int i = 0; i < logdim_; i++)
		{
			rest = dim_ / (offset*divisor);
			B = subfft(B, offset, rest, trafo_);
			offset *= divisor;
		}
		A = B;
		*/
		A.print();
		Reverse(A);
		A.print();
		Reverse(A);
		A.print();
		getchar();
	}

	void ifft(Vectorcd& A)
	{
		Vectorcd B(A.Dim());
		int divisor = 2;
		int offset = 1;
		int rest = 0;
		for (int i = 0; i < logdim_; i++)
		{
			rest = dim_ / (offset*divisor);
			B = subfft(B, offset, rest, trafoback_);
			offset *= divisor;
		}
	}



	Vectorcd subfft(const Vectorcd& A, int offset, int parts, const Vectorcd& twiddle)
	{
		// offset= 2^stage
		// parts = 2^N-stage-1
		Vectorcd B(A.Dim());
		for (int i = 0; i < parts; i++)
		{
			for (int j = 0; j < offset; j++)
			{
				int left   = offset * 2 * j + i;
				int right  = left + offset;
				int twiidx = dim_ / offset*j;
				B(left) = A(left) + twiddle(twiidx)*A(right);
			}
		}
		return B;
	}

protected:
	int dim_;
	int logdim_;
	Vectorcd trafo_;
	Vectorcd trafoback_;
};

