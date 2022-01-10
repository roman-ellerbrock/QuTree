//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorTreeOperator.h"

template<typename T>
Matrix<T> toMatrix(const MLO<T>& M, const Leaf& leaf) {
	Matrix<T> sigma = identityMatrix<T>(leaf.dim());
	for (size_t k = 0; k < M.size(); ++k) {
		if (M.isActive(k, leaf.mode())) {
			const auto& L = M[k];
			const Matrix<T> mk = toMatrix(*L, leaf);
			sigma = mk * sigma;
		}
	}
	return sigma;
}

template <typename T>
void toTensor(Tensor<T>& A, const Matrix<T>& M, size_t part, const Leaf& leaf) {
	for (size_t i = 0; i < A.shape().lastBefore(); ++i) {
		A(i, part) = M[i];
	}
}

template<typename T>
void toTensor(Tensor<T>& A, const MLO<T>& M, size_t part, const Leaf& leaf) {
	/**
	 * Rationale:
	 * Build Tensor representation of the LeafOperator that acts on 'leaf' and store it in A(:,part).
	 * Note: currently only works if there is a single operator acting on 'leaf' in M
	 */
	const TensorShape& shape = A.shape();
	auto sigma = toMatrix(M, leaf);
	toTensor(A, sigma, part, leaf);
}

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& B, size_t l, const Leaf& leaf) {
	size_t dim = leaf.dim();
	Matrix<T> h(dim, dim);
	for (size_t I = 0; I < dim*dim; ++I) {
		h[I] = B(I, l);
	}
	return h;
}

template<typename T>
TensorTreeOperator<T> product(const MLO<T>& M, TensorTreeOperator<T> A, const Tree& tree) {
	for (const Node& node : tree) {
		if(node.isBottomlayer()) {
			const Leaf& leaf = node.getLeaf();
			size_t mode = leaf.mode();
			if (!M.isActive(mode)) { continue; }
			Tensor<T>& phi = A[node];
			const TensorShape& shape = phi.shape();
			///
			Matrix<T> mat = toMatrix(M, leaf);
			///
			for (size_t l = 0; l < shape.lastDimension(); ++l) {
				Matrix<T> hphi = toMatrix(phi, l, leaf);
				hphi = mat * hphi;
				toTensor(phi, hphi, l, leaf);
			}
		}
	}
	return A;
}

template<typename T>
TensorTreeOperator<T>::TensorTreeOperator(const MLO<T>& M,
	const Tree& tree)
	: TensorTreeOperator(tree) {
	occupy(tree);
	vector<size_t> idxs(tree.nLeaves(), 0);

	for (size_t k = 0; k < M.size(); ++k) {
		size_t mode = M.mode(k);
		const Leaf& leaf = tree.getLeaf(mode);
		const auto& node = (const Node&) leaf.parent();

		const shared_ptr<LeafOperator<T>>& h = M[k];
		Matrix<T> hmat = toMatrix(*h, leaf);

		setLeafOperator(hmat, idxs[mode], node);
		idxs[mode]++;
	}
}

template<typename T>
TensorTreeOperator<T>::TensorTreeOperator(const SOP<T>& S,
	const Tree& tree)
	: TensorTreeOperator(tree) {
	occupy(tree);
	assert(S.size() > 0);
}

template<typename T>
TensorTreeOperator<T>::TensorTreeOperator(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorShape& shape = node.shape();
		attributes_.emplace_back(Tensor<T>(shape));
	}
	occupy(tree);
}

template<typename T>
void TensorTreeOperator<T>::occupy(const Tree& tree) {
	for (const Node& node : tree) {
		Tensor<T>& Phi = operator[](node);
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
			setLeafOperator(identityMatrix<T>(dim), 0, node);
		}
	}
}

template<typename T>
void TensorTreeOperator<T>::print(const Tree& tree) const {
	for (const Node& node : tree) {
		node.info();
		attributes_[node.address()].print();
	}
}

template<typename T>
void TensorTreeOperator<T>::setLeafOperator(const Matrix<T>& m,
	size_t operator_idx, const Node& node) {

	const TensorShape& shape = node.shape();
	assert(m.dim1() * m.dim2() == shape.lastBefore());
	assert(operator_idx < shape.lastDimension());

	Tensor<T>& h = operator[](node);
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		h(i, operator_idx) = m[i];
	}
}

template<typename T>
void TensorTreeOperator<T>::occupy(const Tree& tree, mt19937& gen) {
	uniform_real_distribution<double> dist(-1., 1.);
	for (const Node& node : tree) {
		Tensor<T>& A = (*this)[node];
		const TensorShape& shape = A.shape();
		for (size_t i = 0; i < shape.totalDimension(); ++i) {
			A[i] = dist(gen);
		}
		A = qr(A, shape.lastIdx());
	}
}


template class TensorTreeOperator<double>;
template class TensorTreeOperator<complex<double>>;

typedef double d;
typedef complex<double> cd;

template void toTensor(Tensor<d>& A, const MLO<d>& M, size_t part, const Leaf& leaf);
template void toTensor(Tensord& A, const Matrixd& M, size_t part, const Leaf& leaf);
template Matrixd toMatrix(const MLOd& M, const Leaf& leaf);
template Matrixd toMatrix(const Tensord& B, size_t l, const Leaf& leaf);
template TensorTreeOperator<d> product(const MLO<d>& M, TensorTreeOperator<d> A, const Tree& tree);

template void toTensor(Tensorcd& A, const MLOcd& M, size_t part, const Leaf& leaf);
template void toTensor(Tensorcd& A, const Matrixcd& M, size_t part, const Leaf& leaf);
template Matrixcd toMatrix(const MLOcd& M, const Leaf& leaf);
template Matrixcd toMatrix(const Tensorcd& B, size_t l, const Leaf& leaf);
template TensorTreeOperator<cd> product(const MLO<cd>& M, TensorTreeOperator<cd> A, const Tree& tree);

