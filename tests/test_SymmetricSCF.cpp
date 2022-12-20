//
// Created by Roman Ellerbrock on 12/20/22.
//

#include "TreeClasses/SymmetricSCF.h"
#include <gtest/gtest.h>
#include "TreeShape/TreeFactory.h"


TEST (Configuration, Create) {
	Configuration<> c({0, 1});
	Configuration<> d({1, 0});
	ConfigurationTensor<> x({c, d});
	ConfigurationTensor<> y({{0, 1}, {1, 0}});
	ASSERT_EQ(x, y);
}

TEST (Configuration, Product) {
	Configuration<> c({0, 1});
	Configuration<> d({1, 0});
	auto e = c * d;
	ASSERT_EQ(e, Configuration<>({0, 1, 1, 0}));
}

TEST (ConfigurationTensor, Product) {
	ConfigurationTensor<> A({{0, 1}, {1, 0}});
	ConfigurationTensor<> B({{0, 0}, {1, 1}});
	auto C = A * B;
	ASSERT_EQ(C, ConfigurationTensor<>({{0, 1, 0, 0}, {0, 1, 1, 1},
											 {1, 0, 0, 0}, {1, 0, 1, 1}}));
}

TEST (ConfigurationTensor, Sum) {
	ConfigurationTensor<> A({{0, 1}, {1, 0}});
	ConfigurationTensor<> B({{0, 0}, {1, 1}});
	auto C = A + B;
	ASSERT_EQ(C, ConfigurationTensor<>({{0, 1}, {1, 0},
											 {0, 0}, {1, 1}}));
}

TEST (ConfigurationTree, Create) {
	Tree tree = TreeFactory::balancedTree(8, 2, 3);
	ConfigurationTree<> Psi(tree);

	cout << Psi;
}
