//
// Created by Roman Ellerbrock on 8/7/21.
//
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeClasses/TensorTreeFunctions.h"

TensorOperatorTree contractSOP(TensorOperatorTree A, const SOPd& S, const Tree& optree) {

	size_t maxIter = 4;
	double eps = 1e-10;
	cout << "Initial error " << ": ";
	double err = error(A, S, optree);
	cout << err << endl;
	for (size_t i = 0; i < maxIter; ++i) {
		cout << "Iteration: " << i << endl;
		iterate(A, S, optree);
		cout << "Error after iteration " << i << ": ";
		err = error(A, S, optree);
		cout << err << endl;
		if (err < 1e-10) { break; }
	}
	return A;
}

Tensord applyLayer(const TTNOMatrixTree& rep, const TTNOHoleTree& hole,
	const SOPd& S, const Node& node) {

	const TensorShape& shape = node.shape();
	Tensord Bnew(shape);

	if (node.isBottomlayer()) {
		Matrixd shole = hole[node];
		for (size_t l = 0; l < S.size(); ++l) {
			Tensord sterm = S.coeff(l) * toTensor(S, node.getLeaf());
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				auto idx = indexMapping(I, shape);
				size_t i0 = idx[node.parentIdx()];
				Bnew(I) += sterm(I) * shole(i0, l);
			}
		}
	} else {

		vector<Matrixd> Mk = rep.gatherMk(node);
		hole.gatherMk(Mk, node);
/*		for (int i = 0; i < Mk.size(); ++i) {
			cout << "i = " << endl;
			Mk[i].print();
		}*/

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

void iterate(TensorOperatorTree& A, const SOPd& S, const Tree& optree) {
	TTNOMatrixTree rep(S, optree);
	TTNOHoleTree hole(S, optree);
	rep.represent(A, S, optree);
	hole.represent(A, rep, optree);

	for (const Node& node : optree) {
		const Tensord& B = A[node];

		////////////
		if (node.isBottomlayer()){continue;}
//		if (node.isToplayer()){continue;}

		A[node] = applyLayer(rep, hole, S, node);
//		A[node].print();
//		getchar();

		if (!node.isToplayer()) {
			A[node] = qr(A[node]);
			rep.representLayer(A, S, node);
		} else {
			gramSchmidt(A[node]);
			A[node] *= sqrt((double) pow(2, optree.nLeaves()) * S.size()) ;
		}
	}
}

Tensord buildOperator(const MLOd& M, const Leaf& leaf) {

	size_t mode = leaf.mode();
	Matrixd I = identityMatrixd(leaf.dim());
	for (size_t i = 0; i < M.size(); ++i) {
		if (!M.isActive(i, mode)) { continue; }
		const auto& op_ptr = M[i];
		const LeafOperatord& op = *op_ptr;
		Matrixd op_rep = toMatrix(op, leaf);
		I = op_rep * I;
	}
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
		Tensord ml = buildOperator(Ml, leaf);
		Tensord mm = buildOperator(Mm, leaf);
		Matrixd prod = ml.dotProduct(mm);
		x *= abs(prod(0, 0));
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

double error(const TensorOperatorTree& A, const SOPd& S, const Tree& optree) {
	auto AA = TreeFunctions::dotProduct(A, A, optree);
	TTNOMatrixTree AS(S, optree);
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
	err /= normS;
	cout << err << " = " << normA/normS << " - " << overlap/normS << " + " << normS/normS << endl;

	return err;
}

