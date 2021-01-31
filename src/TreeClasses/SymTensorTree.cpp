//
// Created by Roman Ellerbrock on 1/27/21.
//

#include "TreeClasses/SymTensorTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"

void SymTensorTree::initializeFromTT(const TensorTreecd& Psi, const Tree& tree) {
	/// First set up_ to conventional representation (SPFs orthonormal, Toplayer density-weighted)
	up_ = Psi;
	QROrthogonal(up_, tree);
	/// @TODO: Perform simple Dotproduct-test to verify.

	/// Now create orthonormal down-tensors and simultaneously create weighted tensors
	weighted_ = up_;
	down_ = up_;
	QROrthogonalDown(weighted_, down_, up_, tree);
	/// @TODO: Perform simple Contraction-test to verify.

	if (!isWorking(*this, tree)) { cerr << "Normalization of SymTensorTree failed.\n"; exit(1); }


//	rebuild(tree);
}

void SymTensorTree::rebuild(const Tree& tree) {
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		Tensorcd Q = QR(weighted_[node]);
		auto s = Q.DotProduct(up_[node]);
		node.info();
		s.print();
	}
	getchar();
}

Tensorcd solveSLE(const Matrixcd& B, const Tensorcd& A, size_t idx) {
	Matrixcd Aflat = toMatrix(A, idx).Transpose();
	Matrixcd Bt = B.Transpose();
	using namespace Eigen;
	MatrixXcd Bm = Map<MatrixXcd>((complex<double> *) &Bt(0, 0), Bt.Dim1(), Bt.Dim2());
	MatrixXcd Aflatm = Map<MatrixXcd>((complex<double> *) &Aflat(0, 0), Aflat.Dim1(), Aflat.Dim2());

	MatrixXcd x = Bm.colPivHouseholderQr().solve(Aflatm); // QR solver
//	MatrixXcd x = Bm.ldlt().solve(Aflatm); // Alternative solver
	auto r = (Aflatm - Bm * x).norm();
//	cout << "residual:" << r << endl; // Check accuracy
	if (r > 1e-10) {
		cerr << "Normalization of Tensor failed.\n";
		cerr << "Residual: " << r << endl;
		exit(1);
	}

	Matrixcd Abarflat = toQutree(x);
	Abarflat = Abarflat.Transpose();
	Tensorcd Abar = toTensor(Abarflat, A.shape(), idx);
	return Abar;
}

Tensorcd normalizedTensor(const Tensorcd& weighted, size_t k) {
	Matrixcd rho = Contraction(weighted, weighted, k);
	Matrixcd B = toMatrix(sqrt(Diagonalize(rho)));
	Tensorcd Anorm = solveSLE(B, weighted, k);
/*	auto A2 = MatrixTensor(B, Anorm, k);
	if (Residual(weighted, A2) > 1e-7) {
		cerr << "Failed finding normalized tensor!\n";
		getchar();
		exit(1);
	}*/
	return Anorm;
}

template<typename T>
void QROrthogonalDown(TensorTree<T>& weighted, TensorTree<T>& down, const TensorTreecd& up, const Tree& tree) {
	/** Rationale:
	 * a.) Performs an orthogonalization for Top-Down TensorTree.
	 * b.) Has to be in topdown format, i.e. every Tensor is shifted one layer down.
	 * c.) Meant to be used with MatrixTensorTree.
	 * d.) Not tested and probably not working.
	 */

	for (auto it = tree.rbegin(); it != tree.rend(); ++it) {
		const Node& node = *it;
		if (!node.isToplayer()) {
			const Node& parent = node.parent();

			auto rho = Contraction(weighted[parent], weighted[parent], node.childIdx());
			Matrixcd srho = toMatrix(sqrt(Diagonalize(rho)));
			down[node] = solveSLE(srho, weighted[parent], node.childIdx());

/*
			/// Decompose into unitary (Q) and weighted part (R)
			down[node] = QR(weighted[parent], node.childIdx());
			/// Determine Weight from overlap
			Matrixcd B = down[node].DotProduct(weighted[parent]);
			/// Create weighted coefficient
			weighted[node] = MatrixTensor(B, weighted[node], node.nChildren());
			*/
		}
	}
}

bool isWorking(const SymTensorTree& Psi, const Tree& tree) {
	bool works = true;
	double eps = 1e-7;
	/// Check Bottom-up
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.up_[node];
		auto S = Phi.DotProduct(Phi);
		if (Residual(S, IdentityMatrixcd(S.Dim1())) > eps) { works = false; }
	}

	/// Check top-down
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.down_[node];
		auto S = Contraction(Phi, Phi, node.childIdx());
		if (Residual(S, IdentityMatrixcd(S.Dim1())) > eps) { works = false; }
	}

	return works;
}

namespace TreeFunctions {

	void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Phi,
		const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hPhi = TreeFunctions::ApplyHole(hole, Phi, hchild);

		if (!parent.isToplayer() && hole.Active(parent)) {
			hPhi = TensorMatrix(hPhi, hole[parent], parent.childIdx());
		}
		hole[hchild] = Contraction(Phi, hPhi, hchild.childIdx());
	}

	void symContraction(SparseMatrixTreecd& hole, const TensorTreecd& Psi,
		const SparseMatrixTreecd& mat, const Tree& tree) {

		int sub_topnode = mat.Active().size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = mat.Active().MCTDHNode(n);
			if (!node.isToplayer()) {
				symContractionLocal(hole, Psi[node], mat, node);
			}
		}
	}

	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Psi,
		const MLOcd& M, const Tree& tree) {
		TreeFunctions::Represent(mat.first, M, Psi.up_, tree);
		TreeFunctions::symContraction(mat.second, Psi.down_, mat.first, tree);
	}

	void symRepresent(SymMatrixTrees& mats, const SymTensorTree& Psi,
		const SOPcd& S, const Tree& tree) {
		for (size_t l = 0; l < S.size(); ++l) {
			symRepresent(mats[l], Psi, S[l], tree);
		}
	}

	/// Apply Operators

	Tensorcd symApply(const Tensorcd& Phi,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node) {
		Tensorcd hPhi = TreeFunctions::Apply(mats.first, Phi, M, node);
		return symApplyDown(hPhi, mats.second, node);
	}

	void symApply(Tensorcd& HPhi, const Tensorcd& Phi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Node& node) {
		HPhi.Zero();
		for (size_t l = 0; l < H.size(); ++l) {
			HPhi += H.Coeff(l) * symApply(Phi, hmats[l], H[l], node);
		}
	}

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree) {
		for (const Node& node : tree) {
			symApply(HPsi.weighted_	[node], Psi.weighted_[node], hmats, H, node);
		}
	}

}

