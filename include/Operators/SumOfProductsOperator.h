#pragma once
#include "MultiParticleOperator.h"

template<typename T>
class SumOfProductsOperator
	/**
	 * \class SumOfProductsOperator
	 * \ingroup Operators
	 * \brief This class represents a sum-of-products operator.
	 *
	 * SOPs are sums of MultiParticleOperators. The class provides
	 * basic arithmetics that allows to perform high-level operations
	 * on operators.
	 */
{
public:
	// Constructor
	SumOfProductsOperator() = default;

	// Destructor
	~SumOfProductsOperator() = default;

	explicit SumOfProductsOperator(const TTBasis& basis) {
		Initialize(basis);
	}

	explicit SumOfProductsOperator(const MPO<T>& M, T c = 1.);

	void Initialize(const TTBasis& basis) {
		coeff.clear();
		mpos.clear();
		SpecialInitialize(basis);
	}

	// Get the number of MPOs in the Hamiltonian
	virtual int size() const { return mpos.size(); }

	// Get a product-operator from the Hamiltonian
	virtual const MPO<T>& operator()(size_t part) const {
		assert(part < mpos.size());
		return mpos[part];
	}

	// Get a product-operator from the Hamiltonian
	virtual const MPO<T>& operator[](size_t part) const {
		assert(part < mpos.size());
		return mpos[part];
	}

	virtual MPO<T>& operator()(size_t part) {
		assert(part < mpos.size());
		return mpos[part];
	}

	// append a new summand
	void push_back(const MPO<T>& M, complex<double> coeff_) {
		mpos.push_back(M);
		coeff.push_back(coeff_);
	}

	complex<double> Coeff(size_t i) const {
		assert(i < coeff.size());
		return coeff[i];
	}

//	vector<MPO<T>>::const_iterator begin() const {
	auto begin() const {
		return mpos.begin();
	}

//	vector<MPO<T>>::const_iterator end() const {
	auto end() const {
		return mpos.end();
	}

	//////////////////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////////////////
	// multiply with coefficient
	friend SumOfProductsOperator operator*(T c,
		const SumOfProductsOperator& A);

	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		T c);

	// multiply with Multiparticleoperator
	friend SumOfProductsOperator operator*(const MPO<T>& M,
		const SumOfProductsOperator& A);

	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		const MPO<T>& M);

	// multiply with SoP-Operator
	friend SumOfProductsOperator operator*(const SumOfProductsOperator& A,
		const SumOfProductsOperator& B);

	// add SoP-Operator
	friend SumOfProductsOperator operator+(const SumOfProductsOperator& A,
		const SumOfProductsOperator& B);

protected:
	vector<MPO<T>> mpos;
	vector<complex<double> > coeff;

private:
	virtual void SpecialInitialize(const TTBasis& basis) {
		cerr << "Called SpecialInitialize of SOP-base class." << endl;
	}
};

template <typename T>
using SOP = SumOfProductsOperator<T>;

template <typename T>
SOP<T> multAB(const SOP<T>& A, const SOP<T>& B);

