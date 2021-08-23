//
// Created by Roman Ellerbrock on 8/21/21.
//

#include "contractCircuit.h"
#include "TreeClasses/MatrixTensorTreeFunctions.h"

Tensord apply(Tensord A, const MatrixTreed& rep, const MatrixTreed& hole,
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

double error(const TensorTreeOperator& A, const vector<MatrixTreed>& rep,
	const vector<TensorTreeOperator>& SU,
	const SOPd& S, const Tree& optree) {
	const Node& top = optree.topNode();
}

double iterateCircuit(TensorTreeOperator& A, const TensorTreeOperator& U,
	const SOPd& S, const Tree& optree) {
	vector<TensorTreeOperator> SU;
	for (const MLOd& M : S) {
		SU.emplace_back(product(M, U, optree));
	}

	vector<MatrixTreed> rep;
	vector<MatrixTreed> holes;
	for (const TensorTreeOperator& MU : SU) {
		rep.emplace_back(TreeFunctions::dotProduct(A, MU, optree));
		holes.emplace_back(TreeFunctions::contraction(A, MU, rep.back(), optree));
	}
	for (const Node& node : optree) {
		Tensord& phi = A[node];
		phi.zero();
		for (size_t l = 0; l < S.size(); ++l ) {
			Tensord MA = apply(A[node], rep[l], holes[l], node);
			MA *= S.coeff(l);
			phi += MA;
		}
		gramSchmidt(phi);
	}
}


TensorTreeOperator contractCircuit(const SOPd& S, const TensorTreeOperator& U,
	size_t maxIter, const Tree& optree, ostream *os) {

	mt19937 gen(time(nullptr));
	TensorTreeOperator A(optree, gen);
	//	size_t maxIter = 1;
	double eps = 1e-10;
	if (os) *os << "Initial error " << ": ";
	double err = error(U, S, optree);
	if (os) *os << err << endl;
	for (size_t i = 0; i < maxIter; ++i) {
		if (os) *os << "Iteration: " << i << endl;
		err = iterateCircuit(A, U, S, optree);
		if (os) *os << "Error after iteration " << i << ": ";
//		err = error(A, U, S, optree);
		if (os) *os << err << endl;
		if (err < 1e-10) { break; }
	}
	return A;
}


