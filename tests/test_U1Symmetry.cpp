//
// Created by Roman Ellerbrock on 1/25/23.
//

#include "TreeClasses/Discrete/U1Symmetry.h"
#include <gtest/gtest.h>

TEST(U1, CombineLabels) {
	Labels L = {0, 1};
	Labels R = {0, 3, 4};
	Range r(0, 4);
	Labels LR = combine(L, R, r);
	Labels res = {0, 1, 3, 4};
	ASSERT_EQ(LR, res);
}

TEST(U1, partitions) {
	ilist numbers({0, 1, 2, 3});
	auto p = partitions(numbers, 2, 3);
	list<ilist> res({{0, 3}, {1, 2}, {2, 1}, {3, 0}});
	ASSERT_EQ(p, res);
}

TEST(U1, combinatoricalMap) {
	size_t N = 4;
	size_t k = 2;
	vector<ivec> res({
		{1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 1},
		{0, 1, 1, 0}, {0, 1, 0, 1}, {0, 0, 1, 1}});
	vector<ivec> c;
	for (size_t I = 0; I < binomial(k, N); ++I) {
		ivec idx;
		combinatoricalMap(idx, I, k, N);
		c.push_back(idx);
	}
	ASSERT_EQ(res, c);
}

TEST(U1, combinatoricalMapInv) {
	size_t N = 4;
	size_t k = 2;
	vector<ivec> idxs({
		{1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 1},
		{0, 1, 1, 0}, {0, 1, 0, 1}, {0, 0, 1, 1}});
	ivec Is;
	size_t I = 0;
	for (const auto& idx: idxs) {
		combinatoricalMap(I, idx, k, N);
		Is.push_back(I);
	}
	ivec res({0, 1, 2, 3, 4, 5});
	ASSERT_EQ(res, Is);
}

TEST(U1, combinationToPartition) {
	ivec p({0, 1, 1, 0, 1, 1, 1});
	combinationToPartition(p);
	ivec res({0, 2, 3});
//	cout << p << " | " << res << endl;
	ASSERT_EQ(p, res);

	p = {0};
	combinationToPartition(p);
	res = {0, 0};
	ASSERT_EQ(p, res);

	p = {1};
	combinationToPartition(p);
	res = {1};
	ASSERT_EQ(p, res);

	p = {1, 0, 0, 1, 1, 0};
	combinationToPartition(p);
	res = {1, 0, 2, 0};
	ASSERT_EQ(p, res);

	p = {0, 0, 0, 0};
	combinationToPartition(p);
	res = {0, 0, 0, 0, 0};
	ASSERT_EQ(p, res);

	p = {1, 1, 1, 1, 1};
	combinationToPartition(p);
	res = {5};
	ASSERT_EQ(p, res);
}

TEST(U1, partitionMap) {
	size_t nSummands = 2;
	size_t sum = 3;
	size_t n = nPartitions(sum, nSummands);
	ivec p;
	vector<ivec> partitions;
	for (size_t I = 0; I < n; ++I) {
		partitionMap(p, I, sum, nSummands);
		partitions.push_back(p);
	}

	vector<ivec> res({{3, 0}, {2, 1}, {1, 2}, {0, 3}});
	ASSERT_EQ(partitions, res);
}

TEST(U1, partitionsVec) {
	vector<Labels> numbers({{0, 1, 3}, {0, 2, 3}});
	auto p = partitions(numbers, 3);
	vector<Labels> res({{0, 3}, {1, 2}, {3, 0}});
	ASSERT_EQ(p, res);

	numbers = vector<Labels>({{0, 1, 3}, {0, 2, 3}, {1, 2}});
	p = partitions(numbers, 4);
	res = vector<Labels>({{0, 2, 2}, {0, 3, 1}, {1, 2, 1}, {3, 0, 1}});
	ASSERT_EQ(p, res);
}
