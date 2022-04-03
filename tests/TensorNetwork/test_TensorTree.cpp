//
// Created by Roman Ellerbrock on 12/3/21.
//
#include "UnitTest++/UnitTest++.h"
#include "TensorNetwork/contractions.h"
#include "TensorNetwork/TensorTree.h"
#include "Tree/TreeFactory.h"

SUITE (TensorTree) {

	class Trees {
	public:
		Trees() {
			tree_ = balancedTree(4, 3, 2);
		}

		Tree tree_;
	};

	double eps = 1e-10;

	TEST_FIXTURE (Trees, constructorSize) {
		TensorTreecd Psi(tree_);
			CHECK_EQUAL(7, Psi.nodes_.size());
			CHECK_EQUAL(12, Psi.edges_.size());
	}

	TEST_FIXTURE (Trees, TensorShapes) {
		TensorTreecd Psi(tree_, deltacd);
		for (const Node& node : tree_) {
				CHECK_EQUAL(node.shape_, Psi[node].shape_);
		}

		for (const Edge& edge : tree_.edges()) {
				CHECK_EQUAL(edge.from().shape_, Psi[edge].shape_);
		}
	}

	TEST_FIXTURE (Trees, constructorTensor) {
		TensorTreecd Psi(tree_, deltacd);
		for (const Node& node : tree_) {
				CHECK_CLOSE(0., residual(deltacd(node.shape_), Psi[node]), eps);
		}

		for (const Edge& edge : tree_.edges()) {
			const Tensorcd& phi = Psi[edge];
			auto delta = contraction(phi, phi, edge);
				CHECK_CLOSE(0., residual(delta, identitycd(delta.shape_)), eps);
		}
	}

	TEST_FIXTURE (Trees, matrixTree) {
		auto mat = matrixTreecd(tree_);
		for (const Node& node : tree_) {
				CHECK_EQUAL(true, mat.nodes_.contains(node.address_));
			/// @TODO: define how node tensor should look like. Until now there it isn't defined
		}

		for (const Edge& edge : tree_.edges()) {
				CHECK_EQUAL(true, mat.edges_.contains(edge.address()));
				CHECK_EQUAL(edge.shape(), mat[edge].shape_);
		}
	}

	TEST_FIXTURE (Trees, matrixTreeEdges) {
		auto mat = matrixTreecd(tree_, SubTreeParameters({}, false));
		for (const Node& node : tree_) {
				CHECK_EQUAL(false, mat.nodes_.contains(node.address_));
		}

		for (const Edge& edge : tree_.edges()) {
				CHECK_EQUAL(true, mat.edges_.contains(edge.address()));
		}
	}

	TEST_FIXTURE (Trees, matrixTreeNodes) {
		auto mat = matrixTreecd(tree_, SubTreeParameters({}, true, false));
		for (const Node& node : tree_) {
				CHECK_EQUAL(true, mat.nodes_.contains(node.address_));
		}

		for (const Edge& edge : tree_.edges()) {
				CHECK_EQUAL(false, mat.edges_.contains(edge.address()));
		}
	}

	TEST_FIXTURE (Trees, matrixTreePreEdges) {
		auto mat = matrixTreecd(tree_);
		for (const Edge *edge : mat.edges_) {
				CHECK_EQUAL(edge->shape(), mat[edge].shape_);
		}

		for (const Edge *e : mat.edges_) {
			vector<Edge> pre = mat.preEdges(e);
			for (const Edge pe : pre) {
				CHECK_EQUAL(true, mat.edges_.contains(pe.address()));
				CHECK_EQUAL(e->shape(), mat[e].shape_);
			}
		}
	}

	TEST_FIXTURE (Trees, matrixTreePreEdgesSparse) {
		vector<size_t> idx = {1, 3};
		auto mat = matrixTreecd(tree_, idx);
		for (const Edge *e : mat.edges_) {
			auto pre = mat.preEdges(e);
			for (const Edge &edge : pre) {
					CHECK_EQUAL(true, mat.edges_.contains(edge.address()));
					CHECK_EQUAL(edge.shape(), mat[edge].shape_);
			}
		}
	}
}

