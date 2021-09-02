//
// Created by Roman Ellerbrock on 8/21/21.
//

#include "TreeOperators/TensorOperators/contractCircuit.h"
#include "TreeClasses/MatrixTensorTreeFunctions.h"

template<typename T>
Tensor<T> apply(Tensor<T> A, const MatrixTree<T>& rep, const MatrixTree<T>& hole,
	const Node& node) {
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			A = matrixTensor(rep[child], A, k);
		}
	}
	if (!node.isToplayer()) {
		A = matrixTensor(hole[node], A, node.parentIdx());
	}
	return A;
}

template<typename T>
double iterateCircuit(TensorTreeOperator<T>& A, const TensorTreeOperator<T>& U,
	const SOP<T>& S, const Tree& optree) {
	vector<TensorTreeOperator<T>> SU;
	for (const MLO<T>& M : S) {
		SU.emplace_back(product(M, U, optree));
	}

	vector<MatrixTree<T>> rep;
	vector<MatrixTree<T>> holes;
	for (const TensorTreeOperator<T>& MU : SU) {
		rep.emplace_back(TreeFunctions::dotProduct(A, MU, optree));
		holes.emplace_back(TreeFunctions::contraction(A, MU, rep.back(), optree));
	}

	for (const Node& node : optree) {
		Tensor<T>& phi = A[node];
		phi.zero();
		for (size_t l = 0; l < S.size(); ++l) {
			Tensor<T> mu = apply(SU[l][node], rep[l], holes[l], node);
			mu  *= S.coeff(l);
			phi += mu;
		}
		/// QR decomposition and re-represent layer
		if (!node.isToplayer()) {
			phi = qr(phi);
			for (size_t l = 0; l < S.size(); ++l) {
				TreeFunctions::dotProductLocal(rep[l], A[node], SU[l][node], node);
			}
		}
	}
	return 0.;
}

template<typename T>
double error_circuit(const TensorTreeOperator<T>& A, const SOP<T>& S,
	const TensorTreeOperator<T>& U, const Tree& optree) {

	auto Saa = TreeFunctions::dotProduct(A, A, optree)[optree.topNode()];

	vector<TensorTreeOperator<T>> SU;
	for (const MLO<T>& M : S) {
		SU.emplace_back(product(M, U, optree));
	}

	Matrix<T> Sab(optree.nStates(), optree.nStates());
	for (size_t l = 0; l < S.size(); ++l) {
		Matrix<T> summand = TreeFunctions::dotProduct(A, SU[l], optree)[optree.topNode()];
		Sab += S.coeff(l) * summand;
	}

	Matrix<T> Sbb(optree.nStates(), optree.nStates());
	for (size_t i = 0; i < SU.size(); ++i) {
		for (size_t j = 0; j < SU.size(); ++j) {
			Sbb += S.coeff(i)*S.coeff(j)*TreeFunctions::dotProduct(SU[i], SU[j], optree)[optree.topNode()];
		}
	}
	double err = abs(Saa(0, 0)) - 2 * real(Sab(0, 0)) + abs(Sbb(0, 0));
	return err;
}

/**
 * \brief Calculate product of SOP operator with TTNO, i.e. A = S * U
 * @tparam T base-type
 * @param S Sum-of-product operator that is applied
 * @param U contrcted Operator
 * @param maxIter code will leave after max-Iterations if not converged
 * @param optree operator tree
 * @param os out-stream for output during optimization (optional)
 * @return contracted TTNO A
 */
template<typename T>
[[nodiscard]] TensorTreeOperator<T> contractCircuit(const SOP<T>& S, const TensorTreeOperator<T>& U,
	size_t maxIter, const Tree& optree, ostream *os) {

	mt19937 gen(time(nullptr));
	TensorTreeOperator<T> A(optree, gen);
	double eps = 1e-10;
	if (os) *os << "Initial error " << ": ";
	double err = error_circuit(A, S, U, optree);
	if (os) *os << err << endl;
	for (size_t i = 0; i < maxIter; ++i) {
		if (os) *os << "Iteration: " << i << endl;
		iterateCircuit(A, U, S, optree);
		if (os) *os << "Error after iteration " << i << ": ";
		err = error_circuit(A, S, U, optree);
		if (os) *os << err << endl;
		if (err < eps) { break; }
	}
	return A;
}

template<typename T>
TensorTreeOperator<T> contractCircuit(const SOPVector<T>& S, size_t maxIter,
	const Tree& optree, ostream *os) {
	if (S.empty()) {
		return TensorTreeOperator<T>();
	}

	mt19937 gen(time(nullptr));
	TensorTreeOperator<T> A(optree, gen);
	contractSOP(A, S.front(), maxIter, optree, os);

	for (size_t l = 1; l < S.size(); ++l) {
		A = contractCircuit(S[l], A, maxIter, optree, os);
	}
	return A;
}

typedef complex<double> cd;
typedef double d;

template double error_circuit(const TensorTreeOperator<d>& A, const SOP<d>& S,
	const TensorTreeOperator<d>& U, const Tree& optree);
template double error_circuit(const TensorTreeOperator<cd>& A, const SOP<cd>& S,
	const TensorTreeOperator<cd>& U, const Tree& optree);

template double iterateCircuit(TensorTreeOperator<d>& A, const TensorTreeOperator<d>& U,
	const SOP<d>& S, const Tree& optree);
template double iterateCircuit(TensorTreeOperator<cd>& A, const TensorTreeOperator<cd>& U,
	const SOP<cd>& S, const Tree& optree);

template TensorTreeOperator<d> contractCircuit(const SOP<d>& S, const TensorTreeOperator<d>& U,
	size_t maxIter, const Tree& optree, ostream *os);
template TensorTreeOperator<cd> contractCircuit(const SOP<cd>& S, const TensorTreeOperator<cd>& U,
	size_t maxIter, const Tree& optree, ostream *os);

template TensorTreeOperator<d> contractCircuit(const SOPVector<d>& S,
	size_t maxIter, const Tree& optree, ostream *os);
template TensorTreeOperator<cd> contractCircuit(const SOPVector<cd>& S,
	size_t maxIter, const Tree& optree, ostream *os);
