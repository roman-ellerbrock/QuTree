//
// Created by Roman Ellerbrock on 8/7/21.
//
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeClasses/TensorTreeFunctions.h"

TensorTreeOperator contractSOP(TensorTreeOperator A, const SOPd& S,
	size_t maxIter, const Tree& optree, ostream *os) {

//	size_t maxIter = 1;
	double eps = 1e-10;
	if (os) *os << "Initial error " << ": ";
	double err = error(A, S, optree);
	if (os) *os << err << endl;
	for (size_t i = 0; i < maxIter; ++i) {
		if (os) *os << "Iteration: " << i << endl;
		iterate(A, S, optree);
		if (os) *os << "Error after iteration " << i << ": ";
		err = error(A, S, optree);
		if (os) *os << err << endl;
		if (err < 1e-10) { break; }
	}
	return A;
}

Tensord applyLayer(const TTOMatrixTree& rep, const TTOHoleTree& hole,
	const SOPd& S, const Node& node) {

	const TensorShape& shape = node.shape();
	Tensord Bnew(shape);

	if (node.isBottomlayer()) {
		const Matrixd& shole = hole[node];

		for (size_t l = 0; l < S.size(); ++l) {
			Tensord sterm = S.coeff(l) * toTensor(S, node.getLeaf());
			const TensorShape termShape = sterm.shape();
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

		vector<Matrixd> Mk = rep.gatherMk(node);
		hole.gatherMk(Mk, node);
		for (size_t l = 0; l < S.size(); ++l) {
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				auto idx = indexMapping(I, shape);
				double factor = prodMk(idx, Mk, l);
				Bnew(I) += S.coeff(l) * factor;
			}
		}
	}
	return Bnew;
}

void iterate(TensorTreeOperator& A, const SOPd& S, const Tree& optree) {
	TTOMatrixTree rep(S, optree);
	TTOHoleTree hole(S, optree);
	rep.represent(A, S, optree);
	hole.represent(A, rep, optree);

	for (const Node& node : optree) {
		const Tensord& B = A[node];

		////////////
//		if (node.isBottomlayer()) { continue; }
//		if (node.isToplayer()) { continue; }

		A[node] = applyLayer(rep, hole, S, node);

		if (!node.isToplayer()) {
			A[node] = qr(A[node]);
			rep.representLayer(A, S, node);
		} else {
//			A[node] = qr(A[node]);
//			gramSchmidt(A[node]);
//			A[node] *= sqrt((double) pow(2, optree.nLeaves()) * S.size()/2.) ;
		}
	}
}

Tensord buildOperator(const MLOd& M, const Leaf& leaf, bool adjoint) {

	size_t mode = leaf.mode();
	Matrixd I = identityMatrixd(leaf.dim());
	for (size_t i = 0; i < M.size(); ++i) {
		if (!M.isActive(i, mode)) { continue; }
		const auto& op_ptr = M[i];
		const LeafOperatord& op = *op_ptr;
		Matrixd op_rep = toMatrix(op, leaf);
		I = op_rep * I;
	}
	if (adjoint) { I = I.adjoint(); }
	Tensord Itens({leaf.dim() * leaf.dim(), 1});
	for (size_t i = 0; i < Itens.shape().totalDimension(); ++i) {
		Itens[i] = I[i];
	}
	return Itens;
}

double prodnorm(const MLOd& Ml, const MLOd& Mm, const Tree& tree) {
	double x = 1.;
	for (size_t k = 0; k < tree.nLeaves(); ++k) {
		const Leaf& leaf = tree.getLeaf(k);
		Tensord ml = buildOperator(Ml, leaf, false);
		Tensord mm = buildOperator(Mm, leaf, false);
/*		cout << "leaf: " << k << endl;
		ml.print();
		mm.print();*/
		Matrixd prod = ml.dotProduct(mm);
		x *= prod(0, 0);
	}
	return x;
}

double norm(const SOPd& S, const Tree& tree) {
	double norm = 0.;
	for (size_t m = 0; m < S.size(); ++m) {
		const MLOd& Mm = S[m];
		for (size_t l = 0; l < S.size(); ++l) {
			const MLOd& Ml = S[l];
			double factor = prodnorm(Ml, Mm, tree);
			norm += S.coeff(l) * S.coeff(m) * factor;
		}
	}
	return norm;
}

double error(const TensorTreeOperator& A, const SOPd& S, const Tree& optree) {
	auto AA = TreeFunctions::dotProduct(A, A, optree);
	TTOMatrixTree AS(S, optree);
	AS.represent(A, S, optree);

	auto aa = AA[optree.topNode()];
	auto as = AS[optree.topNode()];

	double normA = 0;
	for (size_t j = 0; j < aa.dim2(); ++j) {
		for (size_t i = 0; i < aa.dim1(); ++i) {
			normA += aa(i, j);
		}
	}

	double overlap = 0;
	for (size_t j = 0; j < as.dim2(); ++j) {
		for (size_t i = 0; i < as.dim1(); ++i) {
			overlap += 2 * S.coeff(j) * as(i, j);
		}
	}

	double normS = norm(S, optree);
	double err = normA + normS - overlap;
	if (normS < 1e-15) {
		cerr << "norm of SOP operator too small.\n";
		exit(1);
	}
//	double err = normA + 1 - overlap/sqrt(normS);
//	os << "normS: " << normS << endl;
	err /= normS;
//	os << err << " = " << normA/normS << " - " << overlap / normS << " + " << normS / normS << endl;

	return err;
}

