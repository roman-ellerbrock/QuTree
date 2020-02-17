#include "Util/FFT.h"

size_t FFT::getGoodSize(size_t size)
{
	//check if size is a power of 2
	size_t mod = size % 2;
	if(mod == 0) return size;
	
	//get the next bigger power of 2
	size_t potens = 0;
	while(size > 0)
	{
		size = size/2;
		potens++;
	}
	return pow(2,potens);
}

Tensorcd FFT::reverseOrder(const Tensorcd& in)
{
	//get the next best size for fft
	const TensorDim& dim = in.Dim();
	size_t size = dim.LastBefore();
	size_t states = dim.GetNumTensor();
	size_t N = getGoodSize(size);

	//build the output tensor withe the new sizes
	vector<size_t> d;
	d.push_back(N);
	TensorDim newdim(d, dim.GetNumTensor());
	Tensorcd out(newdim);
	
	for(size_t n = 0; n < states; n++)
		for(size_t i = 0; i < size; i++)
			out[n*size + i] = in[n*size + i];

	complex<double> tmp;

	//reverse the order
	for(size_t n = 0; n < states; n++)
	{
		size_t j = 1;
		for(size_t i = 1; i <= size; i++)
		{
			if(j > i)
			{
				tmp = out[n*size + j-1];
				out[n*size + j-1] = out[n*size + i-1];
				out[n*size + i-1] = tmp;
			}

			size_t m = size/2;
			while((m >= 2) && (j > m))
			{
				j -= m;
				m /= 2;
			}
			j += m;
		} 
	}

	return out;
}

Tensorcd FFT::generalFFT(const Tensorcd& in, const int sign)
{
	//get sizes for fft
	const TensorDim& dim = in.Dim();
	size_t size = dim.LastBefore();
	size_t states = dim.GetNumTensor();

	//get primefactors
	vector<size_t> primefactors = primeFactorization(size);
	size_t numberOfFactors = primefactors.size();

	//make a factor 2 fft if possible
	if(primefactors[numberOfFactors-1] == 2)
		return factor2FFT(in, sign);

	//otherwise perform dft
	return dft(in, sign);
}

Tensorcd FFT::dft(const Tensorcd& in, const int sign)
{
	//get sizes for dft
	const TensorDim& dim = in.Dim();
	size_t size = dim.LastBefore();
	size_t states = dim.GetNumTensor();

	//otherwise perform dft
	vector<size_t> d;
	d.push_back(size);
	TensorDim newdim(d, states);
	Tensorcd out(newdim);
	for(size_t k = 0; k < states; k++)
	{
		for(size_t i = 0; i < size; i++)
		{
			double angle = (2.*M_PI*i*sign)/(1.*size);
			complex<double> factor(cos(angle),sin(angle));
			complex<double> w(1/sqrt(1.*size),0.);
			for(size_t j = 0; j < size; j++)
			{
				out[k*size + i] += w*in[k*size + j];
				w *= factor;
			}
		}
	}
	return out;
}

Tensorcd FFT::factor2FFT(const Tensorcd& in, const int sign)
{
	Tensorcd out = reverseOrder(in);
	
	danielsonLanczosAlgorithm(out, sign);

	return out;
}

void FFT::danielsonLanczosAlgorithm(Tensorcd& in, const int sign)
{
	//get sizes for fft
	const TensorDim& dim = in.Dim();
	size_t size = dim.LastBefore();
	size_t states = dim.GetNumTensor();

	//save some intermediat results
	complex<double> tmp(0., 0.);

	//Transform each vector seperadly
	for(size_t n = 0; n < states; n++)
	{
		//Recursion
		size_t mmax = 1;


		while(mmax < size)
		{
			size_t j = 0;
			size_t step = 2*mmax;
			
			//prefector for exponent
			double angle = M_PI/(1.*sign*mmax);
			complex<double> factor(-2.*pow(sin(0.5*angle),2), sin(angle));

			//transformation factor
			complex<double> w(1., 0.);
			
			//go throug all partitions of the vector
			for(size_t m = 0; m < mmax; m++)
			{
				for(size_t i = m; i < size; i += step)
				{
					//Davidson-Lanczos-Formula
					j = i + mmax;

					tmp = w*in[n*size + j];
				  
					in[n*size + j] = in[n*size + i] - tmp;
					in[n*size + i] = in[n*size + i] + tmp;

				}
				w = w*factor + w;
			}

			mmax = step;
		}
	}

	in *= 1./sqrt(size);
}

vector<size_t> FFT::primeFactorization(size_t number) const
{
	vector<size_t> primefactors;

	if(number == 1)
		primefactors.push_back(1);

	size_t inter = number;
	for(size_t i = 2; i <= number; i++)
	{
		while(inter % i == 0)
		{
			inter /= i;
			primefactors.push_back(i);
		}

		if(inter == 1) break;
	}

	return primefactors;
}
