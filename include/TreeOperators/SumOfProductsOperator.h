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
		initialize(tree);
	}

	explicit SumOfProductsOperator(const MLO<T>& M, T c = 1.);

	void initialize(const Tree& tree) {
		coeff_.clear();
		mpos_.clear();
		specialInitialize(tree);
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

	complex<double> coeff(size_t i) const {
		assert(i < coeff_.size());
		return coeff_[i];
	}

	auto begin() const {
		return mpos_.begin();
	}

	auto end() const {
		return mpos_.end();
	}

	void print(ostream& os = cout) const {
		os << "Number of parts in SOP operator: " << size() << endl;
		for (size_t i = 0; i < size(); ++i) {
			os << "coeff: " << coeff(i) << endl;
			mpos_[i].print(os);
		}
	}

	//////////////////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////////////////
	/// See https://stackoverflow.com/questions/4660123/overloading-friend-operator-for-template-class
	/// for more information
	/// These are extroverts

	template<typename U>
	friend SumOfProductsOperator<U> operator*(U c, const SumOfProductsOperator<U>& A);

	template<typename U>
	friend SumOfProductsOperator<U> operator*(const SumOfProductsOperator<U>& A,
		U c);

	template<typename U>
	friend SumOfProductsOperator<U> operator*(const MLO<U>& M,
		const SumOfProductsOperator<U>& A);

	template<typename U>
	friend SumOfProductsOperator<U> operator*(const SumOfProductsOperator<U>& A,
		const MLO<U>& M);

	template<typename U>
	friend SumOfProductsOperator<U> operator*(const SumOfProductsOperator<U>& A,
		const SumOfProductsOperator<U>& B);

	template<typename U>
	friend SumOfProductsOperator<U> operator+(const SumOfProductsOperator<U>& A,
		const SumOfProductsOperator<U>& B);

protected:
	vector<MLO<T>> mpos_;
	vector<complex<double>> coeff_;

private:
	virtual void specialInitialize(const Tree& tree) {
		cerr << "Called specialInitialize of SOP-base class." << endl;
	}
};

template <typename T>
using SOP = SumOfProductsOperator<T>;

typedef SOP<complex<double>> SOPcd;
typedef SOP<double> SOPd;

template <typename T>
SOP<T> multAB(const SOP<T>& A, const SOP<T>& B);

