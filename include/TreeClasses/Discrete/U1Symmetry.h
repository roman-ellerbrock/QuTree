//
// Created by Roman Ellerbrock on 1/25/23.
//
#ifndef U1SYMMETRY_H
#define U1SYMMETRY_H
#include "stdio.h"
#include <random>
#include "TreeClasses/Discrete/SymmetricSCF.h"

/**
 * Rationale:
 * Reading https://math.stackexchange.com/questions/474741/formula-for-combinations-with-replacement
 * helps understanding the thoughts behind partitions & combinatorics in this file.
 */

size_t factorial(size_t n);

size_t binomial(size_t k, size_t n);

std::vector<size_t> nChooseKrandom(size_t k, size_t n, std::mt19937& gen);

/// Symmetry label, e.g. particle number
using Label = size_t;
using Labels = Configuration<Label>;

/// range of allowed labels, e.g. particle number <= 5 and >= 0
class Range {
public:
	Range(Label min, Label max)
		: min_(min), max_(max) {}

	Label min_{0};
	Label max_{1};

	[[nodiscard]] bool isAllowed(size_t x) const {
		return ((x >= min_) && (x <= max_));
	}

	[[nodiscard]] Labels numbers() const {
		Labels l;
		for (size_t m = min_; m < max_; ++m) {
			l.push_back(m);
		}
		return l;
	}
};

Labels label_combinations(const NodeAttribute<Labels>& up, const NodeAttribute<Labels>& down,
	const Range& range, const Node& node, const Node& hole);

Labels combine(const Labels& L, const Labels& R, const Range& range);

using ilist = list<size_t>;
list<ilist> partitions(ilist numbers, size_t size, size_t sum);
vector<Labels> partitions(const vector<Labels> numbers, size_t sum);
vector<Labels> partitions(const vector<const Labels*> numbers, size_t sum);

using ivec = vector<size_t>;

void combinatoricalMap(ivec& idx, size_t I, size_t k, size_t N);
void combinatoricalMap(size_t& I, const ivec& idx, size_t k, size_t N);
void combinatoricalMapBackconvention(ivec& idx, size_t I, size_t k, size_t N);

void randomCombination(ivec& idx, size_t k, size_t N, mt19937& gen);
ostream& operator<<(ostream& os, const ivec& X);

size_t nPartitions(size_t sum, size_t nSummands);

/// returns the Ith partitions for 'nSummands' summands that sum up to 'sum'.
void partitionMap(ivec& p, size_t& I, size_t sum, size_t nSummands);

/// Helper-function that does: 0111011 -> 032 (interfaced for unit test)
void combinationToPartition(ivec& p);

#endif //U1SYMMETRY_H
