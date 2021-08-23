//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorTreeOperator.h"

Matrixd toMatrix(const MLOd& M, const Leaf& leaf) {
	Matrixd sigma = identityMatrixd(leaf.dim());
	for (size_t k = 0; k < M.size(); ++k) {
		if (M.isActive(k, leaf.mode())) {
			const auto& L = M[k];
			const Matrixd mk = toMatrix(*L, leaf);
			sigma = mk * sigma;
		}
	}
	return sigma;
}

void toTensor(Tensord& A, const Matrixd& M, size_t part, const Leaf& leaf) {
	for (size_t i = 0; i < A.shape().lastBefore(); ++i) {
		A(i, part) = M[i];
	}
}

void toTensor(Tensord& A, const MLOd& M, size_t part, const Leaf& leaf) {
	/**
	 * Rationale:
	 * Build Tensor representation of the LeafOperator that acts on 'leaf' and store it in A(:,part).
	 * Note: currently only works if there is a single operator acting on 'leaf' in M
	 */
	const TensorShape& shape = A.shape();
	auto sigma = toMatrix(M, leaf);
	toTensor(A, sigma, part, leaf);
}

Matrixd toMatrix(const Tensord& B, size_t l, const Leaf& leaf) {
	size_t dim = leaf.dim();
	Matrixd h(dim, dim);
	for (size_t I = 0; I < dim*dim; ++I) {
		h[I] = B(I, l);
	}
	return h;
}

TensorTreeOperator product(const MLOd& M, TensorTreeOperator A, const Tree& tree) {
	for (const Node& node : tree) {
		if(node.isBottomlayer()) {
			const Leaf& leaf = node.getLeaf();
			size_t mode = leaf.mode();
			if (!M.isActive(mode)) { continue; }
			Tensord& phi = A[node];
			const TensorShape& shape = phi.shape();
			///
			Matrixd mat = toMatrix(M, leaf);
			///
			for (size_t l = 0; l < shape.lastDimension(); ++l) {
				Matrixd hphi = toMatrix(phi, l, leaf);
				hphi = mat * hphi;
				toTensor(phi, hphi, l, leaf);
			}
		}
	}
	return A;
}

TensorTreeOperator::TensorTreeOperator(const MLOd& M,
	const Tree& tree)
	: TensorTreeOperator(tree) {
	occupy(tree);
	vector<size_t> idxs(tree.nLeaves(), 0);

	for (size_t k = 0; k < M.size(); ++k) {
		size_t mode = M.mode(k);
		const Leaf& leaf = tree.getLeaf(mode);
		const auto& node = (const Node&) leaf.parent();

		const shared_ptr<LeafOperator<double>>& h = M[k];
		Matrixd hmat = toMatrix(*h, leaf);

		setLeafOperator(hmat, idxs[mode], node);
		idxs[mode]++;
	}
}

TensorTreeOperator::TensorTreeOperator(const SOPd& S,
	const Tree& tree)
	: TensorTreeOperator(tree) {
	occupy(tree);
	assert(S.size() > 0);
}

TensorTreeOperator::TensorTreeOperator(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorShape& shape = node.shape();
		attributes_.emplace_back(Tensord(shape));
	}
	occupy(tree);
}

void TensorTreeOperator::occupy(const Tree& tree) {
	for (const Node& node : tree) {
		Tensord& Phi = operator[](node);
		Phi.zero();
		for (size_t i = 0; i < Phi.shape().lastDimension(); ++i) {
			Phi(i, i) = 1.;
		}
		if (node.isBottomlayer()) {
			if (!(Phi.shape().order() == 2)) {
				cerr << "Wrong order in tensor operator tree.\n";
				exit(1);
			}
			size_t dim = sqrt(1e-10 + Phi.shape().lastBefore());
			setLeafOperator(identityMatrixd(dim), 0, node);
		}
	}
}

void TensorTreeOperator::print(const Tree& tree) const {
	for (const Node& node : tree) {
		node.info();
		attributes_[node.address()].print();
	}
}

void TensorTreeOperator::setLeafOperator(const Matrixd& m,
	size_t operator_idx, const Node& node) {

	const TensorShape& shape = node.shape();
	assert(m.dim1() * m.dim2() == shape.lastBefore());
	assert(operator_idx < shape.lastDimension());

	Tensord& h = operator[](node);
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		h(i, operator_idx) = m[i];
	}
}

void TensorTreeOperator::occupy(const Tree& tree, mt19937& gen) {
	uniform_real_distribution<double> dist(-1., 1.);
	for (const Node& node : tree) {
		Tensord& A = (*this)[node];
		const TensorShape& shape = A.shape();
		for (size_t i = 0; i < shape.totalDimension(); ++i) {
			A[i] = dist(gen);
		}
		A = qr(A, shape.lastIdx());
	}
}

