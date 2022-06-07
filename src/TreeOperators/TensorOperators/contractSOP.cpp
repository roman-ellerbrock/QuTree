//
// Created by Roman Ellerbrock on 8/7/21.
//
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeClasses/TensorTreeFunctions.h"

template<typename T>
void contractSOP(TensorTreeOperator<T>& A, const SOP<T>& S,
	size_t maxIter, const Tree& optree, ostream *os) {

	TTOMatrixTree<T> rep(S, optree);
	TTOHoleTree<T> hole(S, optree);
//	size_t maxIter = 1;
	double eps = 1e-10;
	if (os) *os << "Initial error " << ": ";
//	double err = error(A, S, optree);
	double err = 1.;
	if (os) *os << err << endl;
	for (size_t i = 0; i < maxIter; ++i) {
		if (os) *os << "Iteration: " << i << endl;
		iterate(A, rep, hole, S, optree);
		if (os) *os << "Error after iteration " << i << ": ";
//		err = error(A, S, optree);
		if (os) *os << err << endl;
		if (err < eps) { break; }
	}
}

template<typename T>
Tensor<T> applyLayer(const TTOMatrixTree<T>& rep, const TTOHoleTree<T>& hole,
	const SOP<T>& S, const Node& node) {

	const TensorShape& shape = node.shape();
	Tensor<T> Bnew(shape);

	if (node.isBottomlayer()) {
		const Matrix<T>& shole = hole[node];

		for (size_t l = 0; l < S.size(); ++l) {
			Tensor<T> sterm = S.coeff(l) * toTensor(S, node.getLeaf());
			const TensorShape termShape = sterm.shape();
			auto idx = indexMapping(0, termShape);
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				auto idx = indexMapping(I, termShape);
				size_t i0 = idx[node.parentIdx()];
				/// Evaluate index for sterm
				auto idxS = idx;
				idxS[termShape.lastIdx()] = l;
				size_t L = indexMapping(idxS, termShape);
				Bnew(I) += S.coeff(l) * sterm(L) * shole(i0, l);
			}
		}
	} else {

		vector<Matrix<T>> Mk = rep.gatherMk(node);
		hole.gatherMk(Mk, node);
		auto idx = indexMapping(0, shape);
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			indexMapping(idx, I, shape);
			for (size_t l = 0; l < S.size(); ++l) {
				T factor = prodMk(idx, Mk, l);
				Bnew(I) += S.coeff(l) * factor;
			}
		}
	}
	return Bnew;
}

template<typename T>
void iterate(TensorTreeOperator<T>& A, TTOMatrixTree<T>& rep, TTOHoleTree<T>& hole,
	const SOP<T>& S, const Tree& optree) {

	/// avoid building these before start
	rep.represent(A, S, optree);
	hole.represent(A, rep, optree);

	for (const Node& node : optree) {
		if (node.shape().lastDimension() == node.shape().lastBefore()) { continue; }
		/// Apply layer
		A[node] = applyLayer(rep, hole, S, node);

		///
		if (!node.isToplayer()) {
			/// perform QR-decomposition on Tensor
			A[node] = qr(A[node]);
			/// Normalize bottomlayer with sqrt(dim)
			/// a.) remove norm to build matrices
			if (node.isBottomlayer()) {
				A[node] /= sqrt(sqrt((double) node.shape().lastBefore()));
			}
			/// b.) rebuild representation
			rep.representLayer(A, S, node);
			/// c.) normalize (and account for previous un-normalization)
			if (node.isBottomlayer()) {
				A[node] *= sqrt((double) node.shape().lastBefore());
			}
		}
	}

	/// @TODO: add top-down cycle
}

template<typename T>
Tensor<T> buildOperator(const MLO<T>& M, const Leaf& leaf, bool adjoint) {

	size_t mode = leaf.mode();
	Matrix<T> I = identityMatrix<T>(leaf.dim());
	for (size_t i = 0; i < M.size(); ++i) {
		if (!M.isActive(i, mode)) { continue; }
		const auto& op_ptr = M[i];
		const LeafOperator<T>& op = *op_ptr;
		Matrix<T> op_rep = toMatrix(op, leaf);
		I = op_rep * I;
	}
	if (adjoint) { I = I.adjoint(); }
	Tensor<T> Itens({leaf.dim() * leaf.dim(), 1});
	for (size_t i = 0; i < Itens.shape().totalDimension(); ++i) {
		Itens[i] = I[i];
	}
	return Itens;
}

template<typename T>
T prodnorm(const MLO<T>& Ml, const MLO<T>& Mm, const Tree& tree) {
	T x = 1.;
	for (size_t k = 0; k < tree.nLeaves(); ++k) {
		const Leaf& leaf = tree.getLeaf(k);
		Tensor<T> ml = buildOperator(Ml, leaf, false);
		Tensor<T> mm = buildOperator(Mm, leaf, false);
/*		cout << "leaf: " << k << endl;
		ml.print();
		mm.print();*/
		Matrix<T> prod = ml.dotProduct(mm);
		x *= prod(0, 0);
	}
	return x;
}

template<typename T>
double norm(const SOP<T>& S, const Tree& tree) {
	double norm = 0.;
	for (size_t m = 0; m < S.size(); ++m) {
		const MLO<T>& Mm = S[m];
		for (size_t l = 0; l < S.size(); ++l) {
			const MLO<T>& Ml = S[l];
			T factor = prodnorm(Ml, Mm, tree);
			norm += abs(S.coeff(l) * S.coeff(m) * factor);
		}
	}
	return norm;
}

template<typename T>
double error(const TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& optree, ostream *os) {
	/// @TODO: doesn't work for electronic structure anymore - fix this!

	auto AA = TreeFunctions::dotProduct(A, A, optree);
	TTOMatrixTree<T> AS(S, optree);
	AS.represent(A, S, optree);

	auto aa = AA[optree.topNode()];
	auto as = AS[optree.topNode()];

	double normA = 0;
	for (size_t j = 0; j < aa.dim2(); ++j) {
		for (size_t i = 0; i < aa.dim1(); ++i) {
			normA += abs(aa(i, j));
		}
	}

	double overlap = 0;
	for (size_t j = 0; j < as.dim2(); ++j) {
		for (size_t i = 0; i < as.dim1(); ++i) {
			overlap += real(2. * S.coeff(j) * as(i, j));
		}
	}

//	double err = abs(normA - 0.5 * overlap);
//	return err;

	double normS = norm(S, optree);
//	double normS = 1.;
	double err = normA + normS - overlap;
	if (normS < 1e-15) {
		cerr << "norm of SOP operator too small.\n";
		exit(1);
	}
//	double err = normA + 1 - overlap/sqrt(normS);
//	os << "normS: " << normS << endl;
	err /= normS;
	if (os) { (*os) << err << " = " << normA / normS << " - " << overlap / normS << " + " << normS / normS << endl; }
//	cout << err << " = " << normA << " - " << overlap << " + " << normS << endl;

	return err;
}

typedef double d;

typedef complex<double> cd;

template void contractSOP(TensorTreeOperator<d>& A, const SOP<d>& S,
	size_t maxIter, const Tree& optree, ostream *os);
template void contractSOP<cd>(TensorTreeOperator<cd>& A, const SOP<cd>& S,
	size_t maxIter, const Tree& optree, ostream *os);

template double error(const TensorTreeOperator<d>& A, const SOP<d>& S, const Tree& optree, ostream *os);
template double error(const TensorTreeOperator<cd>& A, const SOP<cd>& S, const Tree& optree, ostream *os);
