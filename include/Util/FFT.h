#pragma once
#include "Core/Tensor.h"

class FFT
{
public:

	FFT() = default;
	~FFT() = default;

	Tensorcd forwardFFT(const Tensorcd& in){return generalFFT(in, -1);}
	Tensorcd backwardFFT(const Tensorcd& in){return generalFFT(in, 1);}
	
protected:
	//fouriertransforms a Tensor statewise
	Tensorcd generalFFT(const Tensorcd& in, const int sign);
	Tensorcd factor2FFT(const Tensorcd& in, const int sign);
	Tensorcd dft(const Tensorcd& in, const int sign);
	
	//rearanges the complex number of each vector
	//The state_ indeces lable the vectors
	Tensorcd reverseOrder(const Tensorcd& in);

	//execute the fft
	void danielsonLanczosAlgorithm(Tensorcd& in, const int sign);

	//calculates the next bigger potens of two
	size_t getGoodSize(size_t size);

	//perform prime factorization
	vector<size_t> primeFactorization(size_t number) const;
};
