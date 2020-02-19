#pragma once
#include "MultiLeafOperator.h"

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

	explicit SumOfProductsOperator(const Tree& tree) {
		Initialize(tree);
	}

	explicit SumOfProductsOperator(const MLO<T>& M, T c = 1.);

	void Initialize(const Tree& tree) {
		coeff_.clear();
		mpos_.clear();
		SpecialInitialize(tree);
	}

	// Get the number of MPOs in the Hamiltonian
	virtual int size() const { return mpos_.size(); }

	// Get a product-operator from the Hamiltonian
	virtual const MLO<T>& operator()(size_t part) const {
		assert(part < mpos_.size());
		return mpos_[part];
	}

	// Get a product-operator from the Hamiltonian
	virtual const MLO<T>& operator[](size_t part) const {
		assert(part < mpos_.size());
		return mpos_[part];
	}

	virtual MLO<T>& operator()(size_t part) {
		assert(part < mpos_.size());
		return mpos_[part];
	}

	// append a new summand
	void push_back(const MLO<T>& M, complex<double> coeff) {
		mpos_.push_back(M);
		coeff_.push_back(coeff);
	}

	complex<double> Coeff(size_t i) const {
		assert(i < coeff_.size());
		return coeff_[i];
	}

	auto begin() const {
		return mpos_.begin();
	}

	auto end() const {
		return mpos_.end();
	}

	//////////////////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////////////////
	// multiply with coefficient
	friend SumOfProductsOperator<T> operator*(T c,
		const SumOfProductsOperator<T>& A);

	friend SumOfProductsOperator<T> operator*(const SumOfProductsOperator<T>& A,
		T c);

	// multiply with Multiparticleoperator
	friend SumOfProductsOperator<T> operator*(const MLO<T>& M,
		const SumOfProductsOperator<T>& A);

	friend SumOfProductsOperator<T> operator*(const SumOfProductsOperator<T>& A,
		const MLO<T>& M);

	// multiply with SoP-Operator
	friend SumOfProductsOperator<T> operator*(const SumOfProductsOperator<T>& A,
		const SumOfProductsOperator<T>& B);

	// add SoP-Operator
	friend SumOfProductsOperator<T> operator+(const SumOfProductsOperator<T>& A,
		const SumOfProductsOperator<T>& B);

protected:
	vector<MLO<T>> mpos_;
	vector<complex<double>> coeff_;

private:
	virtual void SpecialInitialize(const Tree& tree) {
		cerr << "Called SpecialInitialize of SOP-base class." << endl;
	}
};

template <typename T>
using SOP = SumOfProductsOperator<T>;

template <typename T>
SOP<T> multAB(const SOP<T>& A, const SOP<T>& B);

