//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTreeFunctions.h"

namespace TreeFunctions {

	Matrixcd subMatrix(const Tensorcd& B, size_t idx) {
		const TensorShape& shape = B.shape();
		assert(shape.order() == 3);
		size_t dim1 = shape[0];
		size_t dim2 = shape[1];
		assert(idx < shape[2]);
		Matrixcd M(dim1, dim2);
		for (size_t i = 0; i < shape.lastBefore(); ++i) {
			M[i] = B(i, idx);
		}
		return M;
	}

	MatrixList RepresentLower(const Tensorcd& Phi, const Tensorcd& B,
		const MatrixListTree& H, const Node& node) {

		MatrixList hs;
		const TensorShape& Bshape = B.shape();
		const TensorShape& shape = Phi.shape();
		size_t l0 = Bshape.lastDimension();
		for (size_t l = 0; l < l0; ++l) {
			auto mat = subMatrix(B, l);
			auto C = matrixTensor(mat, Phi, 0);
			auto h = contraction(Phi, C, shape.lastIdx());
			hs.emplace_back(h);
		}
		return hs;
	}

	MatrixList RepresentUpper(const Tensorcd& Phi, const Tensorcd& B,
		const MatrixListTree& H, const Node& node) {

		const TensorShape& Bshape = B.shape();
		size_t l0 = Bshape.lastDimension();
		size_t L0 = Bshape.lastBefore();
		const TensorShape& shape = Phi.shape();

		MatrixList hs;
		for (size_t l = 0; l < l0; ++l) {
			Tensorcd C(shape);
			for (size_t L = 0; L < L0; ++L) {
				size_t Ltot = L + L0 * l;
				vector<size_t> Lbreak = indexMapping(Ltot, Bshape);

				Tensorcd Atilde = Phi;
				/// apply hs from below
				for (size_t k = 0; k < node.nChildren(); ++k) {
					const Node& child = node.child(k);
					const MatrixList& subhs = H[child];
					const Matrixcd& subh = subhs[Lbreak[k]];
					Atilde = matrixTensor(subh, Atilde, k);
				}
				Atilde *= B(Ltot);
				C += Atilde;
			}
			hs.emplace_back(contraction(Phi, C, shape.lastIdx()));
		}

		return hs;
	}

	MatrixList Represent(const Tensorcd& Phi, const Tensorcd& B,
		const MatrixListTree& H, const Node& node) {
		vector<Matrixcd> hs;
		if (node.isBottomlayer()) {
			hs = RepresentLower(Phi, B, H, node);
		} else {
			hs = RepresentUpper(Phi, B, H, node);
		}
		return hs;
	}

	MatrixListTree Represent(const TensorTreecd& Psi, const TensorOperatorTree& H,
		const Tree& tree) {

		MatrixListTree Hrep(tree);
		for (const Node& node : tree) {
			Hrep[node] = Represent(Psi[node], H[node], Hrep, node);
		}

		return Hrep;
	}

	MatrixList ContractionUpper(const Tensorcd& Phi, const Tensorcd& B,
		const MatrixListTree& Hrep, const MatrixListTree& Hmean,
		const Node& child) {

		assert(!child.isToplayer());
		const Node& node = child.parent();
		const TensorShape& Bshape = B.shape();
		size_t k = child.childIdx();
		size_t lk = Bshape[k];

		const TensorShape& shape = Phi.shape();
		size_t dim = shape[k];
		MatrixList hmeans(lk, Matrixcd(dim, dim));

		for (size_t L = 0; L < Bshape.totalDimension(); ++L) {
			vector<size_t> Lbreak = indexMapping(L, Bshape);
			Tensorcd Atilde = Phi;
			for (size_t o = 0; o < node.nChildren(); ++o) {
				const Node& otherchild = node.child(o);
				if (o != k) {
					const MatrixList& hs = Hrep[otherchild];
					size_t idx = Lbreak[o];
					Atilde = matrixTensor(hs[idx], Atilde, o);
				}
			}
			if (!node.isToplayer()) {
				const MatrixList& rho = Hmean[node];
				size_t idx = Lbreak.back();
				Atilde = matrixTensor(rho[idx], Atilde, node.parentIdx());
			}
			Atilde *= B[L];
			size_t lkidx = Lbreak[k];
			hmeans[lkidx] += contraction(Phi, Atilde, k);
		}

		return hmeans;
	}

	MatrixListTree Contraction(const TensorTreecd& Psi, const TensorOperatorTree& H,
		MatrixListTree& Hrep, const Tree& tree) {

		MatrixListTree Hmean(tree);
		for (auto it = tree.rbegin(); it != tree.rend(); ++it) {
			const Node& node = *it;
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				Hmean[node] = ContractionUpper(Psi[parent],
					H[parent], Hrep, Hmean, node);
			}
		}

		return Hmean;
	}

}


