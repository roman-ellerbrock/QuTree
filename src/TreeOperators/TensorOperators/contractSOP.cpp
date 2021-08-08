//
// Created by Roman Ellerbrock on 8/7/21.
//
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeClasses/TensorTreeFunctions.h"

TensorOperatorTree contractSOP(const SOPd& S, const Tree& optree) {

	mt19937 gen(time(NULL));
	TensorOperatorTree A(optree, gen);
	size_t maxIter = 10;
	double eps = 1e-10;
	for (size_t i = 0; i < maxIter; ++i) {
		cout << "Iteration: " << i << endl;
		iterate(A, S, optree);
		cout << "Error after iteration " << i << ": ";
		double err = error(A, S, optree);
		cout << err << endl;
		if (err < 1e-10) { break; }
	}
	return A;
}

void iterate(TensorOperatorTree& A, const SOPd& S, const Tree& optree) {
	TTNOMatrixTree rep(S, optree);
	TTNOHoleTree hole(S, optree);
	rep.represent(A, S, optree);
	hole.represent(A, rep, optree);

	qr(A[optree.topNode()]);
	for (const Node& node : optree) {
		const Tensord& B = A[node];
		const TensorShape& shape = B.shape();

		Tensord Bnew(shape);
		vector<Matrixd> Mk = rep.gatherMk(node);
		hole.gatherMk(Mk, node);

		for (size_t l = 0; l < S.size(); ++l) {
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				auto idx = indexMapping(I, shape);
				double factor = prodMk(idx, Mk, l);
				Bnew(I) += S.coeff(l) * B(I) * factor;
			}
		}
		if (!node.isToplayer()) {
			Bnew = qr(Bnew);
			auto s = Bnew.dotProduct(Bnew);
			s.print();
			getchar();

			rep.representLayer(A, S, node);
		}
		A[node] = Bnew;
	}
}

double norm(const SOPd& S, const Tree& tree) {
	double norm = 0.;
//	auto Ssq = S * S; // @TODO: only holds for hermitian matrices
	for (size_t l = 0; l < S.size(); ++l) {
		double c = S.coeff(l);
		const MLOd& M = S[l];
		for (size_t k = 0; k < tree.nLeaves(); ++k) {
			const Leaf& leaf = tree.getLeaf(k);
			Matrixd I = identityMatrixd(sqrt(1e-15 + leaf.dim()));
			if (M.isActive(k, l)) {
				for (size_t i = 0; i < M.size(); ++i) {
					const auto& op_ptr = M[i];
					const auto& op = *op_ptr;
					auto op_rep = toMatrix(op, leaf);
					I = op_rep * I;
				}
			}
			norm += c * I.frobeniusNorm() / (double) leaf.dim();
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

	double err = 0;
	for (size_t j = 0; j < aa.dim2(); ++j) {
		for (size_t i = 0; i < aa.dim1(); ++i) {
			err += aa(i, j);
		}
	}
	err /= (double) pow(2, optree.nLeaves());

	for (size_t j = 0; j < as.dim2(); ++j) {
		for (size_t i = 0; i < as.dim1(); ++i) {
			err -= 2 * S.coeff(j) * as(i, j);
		}
	}

	return 2e-0;
}

