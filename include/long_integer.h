//
// Created by Roman Ellerbrock on 2019-10-22.
//

#ifndef LONG_INTEGER_H
#define LONG_INTEGER_H
#include "stdafx.h"

/**
 * \class long_integer
 * \brief This is a simple class to handle long_integer arithmetic
 *
 * This class is designed to handle modulo-related arithmetic.
 * Its original purpose has been to deal with the classical part of
 * Shor's algorithm. It is therefore not expected to be especially fast
 * but rather to do the operations without falling into some pitfalls.
 * Shor's algorithm features bitshift and modulo combined arithmetic which
 * can be done simultaneously.
 *
 */

class long_integer {
public:
	explicit long_integer(size_t a = 0, size_t n_bit = 32);

	void add(long_integer a);
	void substract(long_integer a);

	void mod(const long_integer& N);

	bool largest() const { return bits.front(); }

	bool smallest() const { return bits.back(); }

	/// Shift this integer by x bits, i.e. multiply with 2**x.
	void bitshift(size_t i);

	/// Perform bitshift and modulo operation (mod is performed for each bitshift independently)
	void bitshift_mod(size_t x, const long_integer& N);

	bool smaller_than(long_integer a) const {
		size_t n_bit = max(a.size(), size());
		a.pad(n_bit);
		long_integer b(*this);
		b.pad(n_bit);
		return b.smaller_than_samesize(a);
	}

	bool smaller_than(size_t a) const {
		return smaller_than_samesize(long_integer(a, this->size()));
	}

	bool smaller_than_samesize(const long_integer& a) const;

	bool equals(const long_integer& a) const;

	bool equals(size_t a) const { return equals(long_integer(a, 32)); }

	void set_bit(size_t i, bool bit) {
		assert(i < bits.size());
		bits[i] = bit;
	}

	bool operator[](size_t i) const {
		return bits[i];
	}

	void pad(size_t n_bit);
	void cromp(size_t n_bit);
	size_t convert() const;

	void print() const;
	void print_4byte() const;

	size_t size() const { return bits.size(); }

private:

	vector<bool> bits;
};

bool carry(bool a, bool b, bool c);
bool sum(bool b, bool a, bool c);

bool eval_borrow(bool b, bool a, bool borrow);
bool subst(bool b, bool a, bool borrow);

long_integer mult(long_integer a, long_integer b);
long_integer mult_mod(const long_integer& a, const long_integer& b, long_integer N);

size_t gcd(size_t a, size_t b);
void gcdx(size_t& g, size_t& x, size_t& y, size_t a, size_t b);
long_integer EuclidicInverse(const long_integer& a, const long_integer& N);

void align(long_integer& a, long_integer& b);

#endif //LONG_INTEGER_H
