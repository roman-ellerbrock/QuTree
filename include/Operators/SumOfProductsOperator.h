#pragma once
#include "MultiParticleOperator.h"

/**
 * \defgroup Operators
 * \brief This group is a bundle of classes that handle
 * operators in mctdh++. System specific operators are NOT included here.
 *
 * This group contains general operator classes that are required to
 * evaluate equations of motion of MCTDH.
 *
 */

class SumOfProductsOperator
{
public:
	// Constructor
	SumOfProductsOperator() = default;

	// Destructor
	~SumOfProductsOperator() = default;

	SumOfProductsOperator(const mctdhBasis& basis) {
		Initialize(basis);
	}

	SumOfProductsOperator(const MultiParticleOperator& M, complex<double> c = 1.) {
	    push_back(M, c);
	}

	void Initialize(const mctdhBasis& basis) {
		coeff.clear();
		mpos.clear();
		SpecialInitialize(basis);
	}

	// Get the number of MPOs in the Hamiltonian
	virtual int size()const { return mpos.size(); }

	// Get a product-operator from the Hamiltonian
	virtual const MultiParticleOperator& operator()(size_t part)const
	{
		assert(part < mpos.size());
		return mpos[part];
	}

	// Get a product-operator from the Hamiltonian
    virtual const MultiParticleOperator& operator[](size_t part)const
    {
        assert(part < mpos.size());
        return mpos[part];
    }

	virtual MultiParticleOperator& operator()(size_t part)
	{
		assert(part < mpos.size());
		return mpos[part];
	}

	// append a new summand
	void push_back(const MultiParticleOperator& M, complex<double> coeff_)
	{
		mpos.push_back(M);
		coeff.push_back(coeff_);
	}

	complex<double> Coeff(size_t i)const
	{
		assert(i < coeff.size());
		return coeff[i];
	}

	vector<MultiParticleOperator>::const_iterator begin() const{
		return mpos.begin();
	}

	vector<MultiParticleOperator>::const_iterator end() const{
		return mpos.end();
	}

	//////////////////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////////////////
	// multiply with coefficient
	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		complex<double> c)
	{
		SumOfProductsOperator C;
		for (size_t i = 0; i < A.size(); i++)
		{
			C.push_back(A(i), A.Coeff(i) * c);
		}
		return C;
	}

	friend SumOfProductsOperator operator*(const complex<double> c,
		const SumOfProductsOperator& A)
	{
		SumOfProductsOperator C;
		for (size_t i = 0; i < A.size(); i++)
		{
			C.push_back(A(i), c * A.Coeff(i));
		}
		return C;
	}

	// multiply with Multiparticleoperator
	friend SumOfProductsOperator operator*(const MultiParticleOperator& M,
		const SumOfProductsOperator& A)
	{
		SumOfProductsOperator C;
		for (size_t i = 0; i < A.size(); i++)
		{
			MultiParticleOperator MA = M * A(i);
			C.push_back(MA, A.Coeff(i));
		}
		return C;
	}

	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		const MultiParticleOperator& M)
	{
		SumOfProductsOperator C;
		for (size_t i = 0; i < A.size(); i++)
		{
			MultiParticleOperator MA = A(i) * M;
			C.push_back(MA, A.Coeff(i));
		}
		return C;
	}

	// multiply with SoP-Operator
	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		const SumOfProductsOperator& B)
	{
		SumOfProductsOperator C;
		for (size_t i = 0; i < A.size(); i++)
		{
			const MultiParticleOperator Ai = A(i);
			auto acoeff = A.Coeff(i);
			for (size_t j = 0; j < B.size(); j++)
			{
				const MultiParticleOperator Bj = B(j);
				auto bcoeff = B.Coeff(j);
				auto ccoeff = acoeff*bcoeff;
				C.push_back(Ai*Bj, ccoeff);
			}
		}
		return C;
	}

	// add SoP-Operator
	friend SumOfProductsOperator operator+(const SumOfProductsOperator& A,
		const SumOfProductsOperator& B)
	{
		SumOfProductsOperator C = B;

		for (size_t i = 0; i < A.size(); i++)
		{
			C.push_back(A(i), A.Coeff(i));
		}

		return C;
	}

protected:
	vector< MultiParticleOperator > mpos;
	vector< complex<double> > coeff;

private:
	virtual void SpecialInitialize(const mctdhBasis& basis) {
		cerr << "Called SpecialInitialize of SOP-base class." << endl;
	}
};

typedef SumOfProductsOperator SOP;

SumOfProductsOperator multAB(const SumOfProductsOperator& A,
	const SumOfProductsOperator& B);

