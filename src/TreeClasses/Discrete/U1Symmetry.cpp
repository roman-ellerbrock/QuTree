//
// Created by Roman Ellerbrock on 1/25/23.
//
#include "TreeClasses/Discrete/U1Symmetry.h"
#include <list>

using namespace std;


size_t factorial(size_t n) {
	size_t res{1};
	for (size_t k = 1; k <= n; ++k) {
		res *= k;
	}
	return res;
}

size_t binomial(size_t k, size_t n) {
	return factorial(n) / (factorial(n - k) * factorial(k));
}

size_t sumList(const ilist& l) {
	size_t s = 0;
	for (auto x: l) {
		s += x;
	}
	return s;
}

size_t sumVec(const Labels& l) {
	size_t s = 0;
	for (auto x: l) {
		s += x;
	}
	return s;
}

void helper(list<ilist>& res, ilist numbers, ilist slate, size_t size, size_t sum) {
	if (slate.size() > size) { return; }
	if (sumList(slate) == sum && (slate.size() == size)) {
		res.push_back(slate);
		return;
	}

	for (auto x: numbers) {
		slate.push_back(x);
		helper(res, numbers, slate, size, sum);
		slate.pop_back();
	}
}

list<ilist> partitions(ilist numbers, size_t size, size_t sum) {
	list<ilist> res;
	ilist slate;
	helper(res, numbers, slate, size, sum);
	return res;
}

void helper(vector<Labels>& res, const vector<Labels>& numbers, size_t depth, Labels slate, size_t size, size_t sum) {
	if (slate.size() > size) { return; }
	if (sumVec(slate) == sum && (slate.size() == size)) {
		res.push_back(slate);
		return;
	}

	if (depth >= numbers.size()) { return; }
	for (auto x: numbers[depth]) {
		slate.push_back(x);
		size_t ndepth = depth + 1;
		helper(res, numbers, ndepth, slate, size, sum);
		slate.pop_back();
	}
}

vector<Labels> partitions(const vector<Labels> numbers, size_t sum) {
		vector<Labels> res;
		Labels slate;
		size_t depth = 0;
		size_t size = numbers.size();
		helper(res, numbers, depth, slate, size, sum);
		return res;
}

vector<Labels> partitions(const vector<const Labels*> numbers, size_t sum) {
	vector<Labels> nums;
	for (auto x : numbers) {
		nums.push_back(*x);
	}
	return partitions(nums, sum);
}

ostream& operator<<(ostream& os, const ivec& X) {
	for (size_t i = 0; i < X.size(); ++i) {
		os << X[i] << ", ";
	}
	return os;
}

void combinatoricalMap(ivec& idx, size_t I, size_t& now, size_t k, size_t N) {
	/// catch trivial case & sanity check
	if (N == 0) { return; }
	size_t s = 0;
	if (k > 0) { s = binomial(k - 1, N - 1); }
	if (I < s) {
		idx[now] = 1;
		combinatoricalMap(idx, I, ++now, k - 1, N - 1);
	} else {
		idx[now] = 0;
		I -= s;
		combinatoricalMap(idx, I, ++now, k, N - 1);
	}
}

void combinatoricalMap(ivec& idx, size_t I, size_t k, size_t N) {
	size_t now = 0;
	idx.resize(N);
	combinatoricalMap(idx, I, now, k, N);
}

void combinatoricalMap(size_t& I, const ivec& idx, size_t& now, size_t k, size_t N) {
	if (now >= idx.size()) { return; }
	if (idx[now]) {
		combinatoricalMap(I, idx, ++now, k - 1, N - 1);
	} else {
		size_t s = 0;
		if (k > 0) { s = binomial(k - 1, N - 1); }
		I += s;
		combinatoricalMap(I, idx, ++now, k, N - 1);
	}
}

void combinatoricalMap(size_t& I, const ivec& idx, size_t k, size_t N) {
	size_t now = 0;
	I = 0;
	combinatoricalMap(I, idx, now, k, N);
}

ivec nChooseKrandom(size_t k, size_t n, mt19937& gen) {
	ivec vec(n);
	std::iota(std::begin(vec), std::end(vec), 0); //0 is the starting number
	std::shuffle(vec.begin(), vec.end(), gen);
	vec.resize(k);
	return vec;
}


void combinationToPartition(ivec& p) {
	/// 0111011 -> 032
	size_t sum = 0;
	size_t n = 1;
	size_t now = 0;
	for (size_t i = 0; i < p.size(); ++i) {
		if (p[i] == 0) {
			p[now++] = sum;
			sum = 0;
			n++;
		} else {
			sum++;
		}
	}
	p[now] = sum;
	p.resize(n);
}

void partitionMap(ivec& idx, size_t& I, size_t sum, size_t nSummands) {
	/// e.g. sum = 3, nSummands = 2
	/// 3 = 0 + 3 = 1 + 2 = 2 + 1 = 3 + 0
	size_t N = sum + nSummands - 1;
	size_t k = sum;
	idx.resize(N);
	combinatoricalMap(idx, I, k, N);
	combinationToPartition(idx);
}

size_t nPartitions(size_t sum, size_t nSummands) {
	size_t N = sum + nSummands - 1;
	size_t k = sum;
	return binomial(k, N);
}
