//
// Created by Roman Ellerbrock on 2019-10-22.
//
#include "long_integer.h"

vector<bool> padded(vector<bool> a, size_t n_bit) {
	if (n_bit > a.size()) {
	}
	assert(a.size() <= n_bit);
	for (size_t i = a.size(); i < n_bit; ++i) {
		a.push_back(0);
	}
	assert(a.size() == n_bit);
	return a;
}

// Returns binary rep of a in big-endian form
vector<bool> binary_rep(size_t a, size_t n_bit) {
	vector<bool> bin_rep;
	while (a > 0) {
		bin_rep.push_back(a % 2);
		a /= 2;
	}
	bin_rep = padded(bin_rep, n_bit);
	return bin_rep;
}

long_integer::long_integer(size_t a, size_t n_bit) {
	bits = binary_rep(a, n_bit);
}

void long_integer::bitshift(size_t x) {
	/// Shift this integer by x bits, i.e. multiply with 2**x.
	assert(x < size());
	for (int i = size() - 1; i >= x; --i) {
		bits[i] = bits[i - x];
	}
	for (int i = 0; i < x; i++) {
		bits[i] = false;
	}
}

void long_integer::bitshift_mod(size_t x, const long_integer& N) {
	for (size_t i = 0; i < x; ++i) {
		bitshift(1);
		mod(N);
	}
}

void align(long_integer& a, long_integer& b) {
	if (a.size() == b.size()) { return; }
	a.pad(b.size());
	b.pad(a.size());
}

void long_integer::cromp(size_t n_bit) {
	vector<bool> newbit;
	int n_ommit = size() - n_bit;
	if (n_ommit <= 0) { return; }
	vector<bool>::const_iterator first = bits.begin();
	vector<bool>::const_iterator last = bits.end() - n_ommit;
	bits = vector<bool>(first, last);
}

void long_integer::pad(size_t n_bit) {
	for (size_t i = size(); i < n_bit; ++i) {
		bits.push_back(false);
	}
}

void long_integer::mod(const long_integer& N) {
	while (!smaller_than_samesize(N)) {
		substract(N);
	}
}

bool long_integer::equals(const long_integer& a) const {
	assert(a.size() == size());
	for (size_t i = 0; i < a.size(); ++i) {
		if (bits[i] != a[i]) {
			return false;
		}
	}
	return true;
}

bool long_integer::smaller_than_samesize(const long_integer& a) const {
	assert(a.size() == size());
	for (int i = a.size() - 1; i >= 0; --i) {
		if (a[i] && !bits[i]) {
			return true;
		}
		if (bits[i] && !a[i]) {
			return false;
		}
	}
	return false;
}

void long_integer::add(long_integer a) {
	align(a, *this);

	bool c = false;
	bool last_carry = false;
	for (size_t i = 0; i < a.size() - 1; ++i) {
		c = carry(a[i], bits[i], last_carry);
		bits[i] = sum(bits[i], a[i], last_carry);
		last_carry = c;
	}
}

void long_integer::substract(long_integer a) {
	// b <- b - a where b = *this
	align(a, *this);

	bool borrow = false;
	bool next_borrow = false;
	for (size_t i = 0; i < a.size() - 1; ++i) {
		next_borrow = eval_borrow(bits[i], a[i], borrow);
		bits[i] = subst(bits[i], a[i], borrow);
		borrow = next_borrow;
	}
}

bool subst(bool b, bool a, bool borrow) {
	if (borrow) {
		b = !b;
	}
	return a ^ b;
}

bool eval_borrow(bool b, bool a, bool borrow) {
	/// If there is a borrow bit but nothing to borrow from, borrow from next
	if ((!b) && borrow) {
		return true;
	}
	/// Otherwise borrow from b;
	if (borrow) {
		b = !b;
	}
	/// If b empty but a full then borrow from next
	if ((!b) && a) {
		return true;
	}
	/// Otherwise we don't need to borrow
	return false;
}

bool carry(bool a, bool b, bool c) {
	bool carry = false;
	if (b && c) {
		carry = !carry;
	}
	if (b) {
		c = !c;
	}
	if (a && c) {
		carry = !carry;
	}
	return carry;
}

bool sum(bool b, bool a, bool c) {
	if (a) {
		b = !b;
	}
	if (c) {
		b = !b;
	}
	return b;
}

void long_integer::print() const {
	cout << "number:\n";
	for (size_t i = 0; i < size(); ++i) {
		cout << bits[size() - 1 - i] << " ";
	}
	cout << endl;
}

size_t long_integer::convert() const {
	size_t a = 0;
	size_t n_bit = min(32, (int) size());
	for (size_t i = 0; i < n_bit; ++i) {
		if (bits[i]) {
			a += (size_t) pow(2, i);
		}
	}
	return a;
}

void long_integer::print_4byte() const {
	cout << "4 byte int representation: " << convert() << endl;
}

size_t gcd(size_t a, size_t b) {
	size_t g = 0;
	size_t x = 0;
	size_t y = 0;
	gcdx(g, x, y, a, b);
	return g;
}

long_integer mult(long_integer a, long_integer b) {
	align(a, b);
	long_integer c(0, 2 * a.size());
	/// How many bits shifted?
	for (size_t i = 0; i < a.size(); ++i) {
		if (b[i]) {
			long_integer t(0, 2 * a.size());
			for (size_t j = 0; j < a.size() - i; ++j) {
				/// shift a by i bit
				t.set_bit(j + i, a[j]);
			}
			c.add(t);
		}
	}
	return c;
}

long_integer mult_mod(const long_integer& a, const long_integer& b, long_integer N) {
	size_t original_nbit = a.size();
	auto m = mult(a, b);
	align(m, N);
	m.mod(N);
	m.cromp(original_nbit);
	return m;
}

void gcdx(size_t& g, size_t& x, size_t& y, size_t a, size_t b) {
	size_t x0 = 0;
	size_t x1 = 1;
	size_t y0 = 1;
	size_t y1 = 0;
	while (a != 0) {
		size_t q = b / a;
		size_t tmp = a;
		a = b % a;
		b = tmp;
		tmp = y0;
		y0 = y1;
		y1 = tmp - q * y1;
		tmp = x0;
		x0 = x1;
		x1 = tmp - q * x1;
	}
	g = b;
	x = x0;
	y = y0;
}

long_integer EuclidicInverse(const long_integer& a, const long_integer& N) {
	assert(a.size() <= 32);
	assert(N.size() <= 32);
	size_t aa = a.convert();
	size_t NN = N.convert();
	size_t g = 0;
	size_t x = 0;
	size_t y = 0;
	gcdx(g, x, y, aa, NN);
	size_t max = 0;
	max -= NN;
	if (x > max) { x += NN; }
	if (y > max) { y += NN; }
	cout << "g = " << g << endl;
	cout << "x = " << x << endl;
	cout << "y = " << y << endl;
	if (g != 1) {
		cerr << "gcd(a, N) is not 1!\nPlease choose valid numbers.\n";
	}
	return long_integer(x, a.size());
}
