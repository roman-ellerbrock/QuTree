//
// Created by Roman Ellerbrock on 3/11/20.
//

#ifndef XMATRIXTREES_H
#define XMATRIXTREES_H
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Core/Tensor_Extension.h"

using Wavefunction = TensorTreecd;

SOPcd Xsop(const Tree& tree);

Wavefunction Regularize(Wavefunction Psi, const Tree& tree, double eps);

typedef pair<Matrixcd, size_t> MatrixIdx;

class XMatrixTrees {
public:
	explicit XMatrixTrees(const Tree& tree)
		: xops_(Xsop(tree)) , eps_(1e-8) {
		LeafFuncd x = &LeafInterface::applyX;

		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.getLeaf(l);
			size_t mode = leaf.mode();
			MLOcd M(x, mode);

			mats_.emplace_back(SparseMatrixTreecd(M, tree));
			holes_.emplace_back(SparseMatrixTreecd(M, tree, false, true));
		}

	}

	~XMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		using namespace TreeFunctions;
		assert(xops_.size() == mats_.size());
		assert(xops_.size() == holes_.size());
		/**
		 * Analysis:
		 * - density is not included in mean-field x-matrices of lower layers
		 *   in ml trees for inverse trees.
		 * Find Evidence:
		 * - Calculate full mean-field x-matrices (dense tree)
		 *
		 */
		represent(mats_, xops_, Psi, Psi, tree);
		Wavefunction Chi = Regularize(Psi, tree, sqrt(eps_));
		auto rho = TreeFunctions::contraction(Chi, tree, true);
		contraction(holes_, Chi, Chi, mats_, rho, tree);
		UnweightContractions(holes_, Chi, tree);
	}

	void UnweightContractions(vector<SparseMatrixTreecd>& holes,
		const Wavefunction& Psi, const Tree& tree) const {
		auto rho = TreeFunctions::contraction(Psi, tree, true);
		auto rho_sqrt = sqrt(rho, tree);
		auto isqrt_rho = inverse(rho_sqrt, tree, eps_);

		for (auto& xhole : holes) {
			const auto& stree = xhole.sparseTree();
			for (const Node *node_ptr : stree) {
				const Node& node = *node_ptr;
				if (!node.isToplayer()) {
					const auto& isq_rho = isqrt_rho[node];
					auto& x = xhole[node];
//					node.info();
//					x.print();
//					isq_rho.print();
					x = isq_rho * xhole[node] * isq_rho;
//					x.print();
				}
			}
		}
	}

	Tensorcd Optimize(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node, const Node& node_small) const;

	Wavefunction Optimize(Wavefunction Psi,
		const MatrixTreecd& rho, const Tree& tree, const Tree& tree_small);

	void print() const {
		cout << "Xs:\n";
		for (const auto& x : mats_) {
			x.print();
		}

		cout << "X holes:\n";
		for (const auto& xhole : holes_) {
			xhole.print();
		}
	}

	size_t size() const {
		size_t n = xops_.size();
		assert(n == mats_.size());
		assert(n == holes_.size());
		return n;
	}

	SOPcd xops_;
	SparseMatrixTreescd mats_;
	SparseMatrixTreescd holes_;

	double eps_;
};

Matrixcd UnProject(size_t n_occupied, const Matrixcd& X,
	const Tensorcd& Phi);

Tensorcd Occupy(const Tensorcd& Phi, const Matrixcd& trafo,
	size_t n_occupied, const Node& node);

Matrixcd BuildX(const Tensorcd& Phi, const Matrixcd& rho,
	const SparseMatrixTreescd& xmats, const Node& node, double eps);

#endif //XMATRIXTREES_H
