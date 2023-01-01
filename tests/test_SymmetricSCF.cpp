//
// Created by Roman Ellerbrock on 12/20/22.
//

#include "TreeClasses/SymmetricSCF.h"
#include <gtest/gtest.h>
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Extension.h"
#include "TreeClasses/TreeIO.h"
#include "Util/QMConstants.h"

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

TEST (ConfigurationTensor, Random) {
	mt19937 gen(0); /// create predictable series
	TensorShape shape({10, 3});
	ConfigurationTensor<> A = randomConfigurationTensor(shape, gen);
	ASSERT_EQ(A, ConfigurationTensor<>({{0}, {1}, {5}}));
}

TEST (ConfigurationTree, Create) {
	mt19937 gen(1); /// create predictable series
	Tree tree = TreeFactory::balancedTree(4, 2, 3);
	auto Psi = randomConfigurationTree(tree, gen);
}

TEST (ConfigurationTensor, botton) {
	TensorShape shape({10, 5});
	auto A = bottomTensor(shape);
}

TEST (Configuration, to_integer) {
	Configuration<> c({1, 0, 1});
	size_t x = to_integer(c);
	ASSERT_EQ(5, x);

	Configuration<> d({0, 0, 1, 1, 1});
	x = to_integer(d);
	ASSERT_EQ(28, x);
}

TEST (Configuration, split_integers) {
	Configuration<> c({1, 0, 1, 0, 1, 0});
	/// 5 & 2
	auto xs = split_integers(c, 2);
	ASSERT_EQ(5, xs[0]);
	ASSERT_EQ(2, xs[1]);
}

TEST (ConfigurationTree, OptimizePolynomial) {
	mt19937 gen(1); /// create predictable series
	size_t N = 10;

	Tree tree = TreeFactory::balancedTree(N, 2, 4, 2);
	auto Psi = randomConfigurationTree(tree, gen);

	auto f = [](const Configuration<>& c) {
		auto xs = split_integers(c, 2);
		double x = xs[0] - 5.;
		double y = xs[1] - 2.;
		double fx = x * x + y * y;
//		cout << xs[0] << " - " << xs[1] << " | " << fx << endl;
		return fx;
	};

	auto x = optimize(Psi, f, tree, 10, 0);
	ASSERT_EQ(0, f(x));
}

TEST (ConfigurationTree, OptimizeSin) {
	mt19937 gen(1); /// create predictable series
	size_t n = 5; /// qubits/integer
	size_t N = 2 * n;

	Tree tree = TreeFactory::balancedTree(N, 2, 8, 2);
	auto Psi = randomConfigurationTree(tree, gen);

	auto f = [n](const Configuration<>& c) {
		auto x = split_doubles(c, 2);
		double fx = sin(x[0] * QM::two_pi) + cos(x[1] * QM::two_pi) + 2.;
		return fx;
	};

	auto c = optimize(Psi, f, tree, 10, 1);
	ASSERT_NEAR(0., f(c), 1e-4);
}

TEST (ConfigurationTree, Optimize) {
	mt19937 gen(1); /// create predictable series
	size_t N = 12;

	Tree tree = TreeFactory::balancedTree(N, 2, 128, 2);
	auto Psi = randomConfigurationTree(tree, gen);

	/**
	 * 2		-115.7
	 * 3		-125.1
	 * 5		-127.9
	 * 7
	 * 10		-223
	 * 20		-118.3		51915951
	 *
	 * 3	6	-215.409
	 * 3	6	-220.887 (integral)
	 */
/*	for (const Node& node: tree) {
		node.info();
		node.shape().print();
	}
	getchar();*/

	Matrixd Q(N, N);
	Matrixd mu(N, N);
	Tensor_Extension::generate(Q, gen);
	Tensor_Extension::generate(mu, gen);
//	Q = 0.5 * (Q + Q.adjoint());
	auto f = [Q, mu](const Configuration<>& c) {
		Tensord y({c.size()});
		for (size_t j = 0; j < c.size(); ++j) {
			for (size_t i = 0; i < c.size(); ++i) {
				y[j] += Q[i + j * c.size()] * c[i];
			}
		}

		double e = 0.;
		for (size_t i = 0; i < c.size(); ++i) {
			e += c[i] * y[i];
		}
/*		Tensord y = matrixTensor(Q, x, 0);
		auto E = contraction(x, y, 1);
		double e = E[0];*/
/*		for (size_t i = 0; i < c.size(); ++i) {
			e -= mu[i] * c[i];
		}*/
		return e;
	};

/*	auto f = [](const Configuration<>& c) {
		Vectord x(c.size());
		double e = 0.;
		for (size_t i = 0; i < c.size(); ++i) {
			e += c[i];
		}
		return e;
	};*/

//	optimize(Psi, f, tree, 10);
//	getchar();
}