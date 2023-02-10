//
// Created by Roman Ellerbrock on 12/2/20.
//

#ifndef DELTAOPERATOR_H
#define DELTAOPERATOR_H
#include "Core/Matrix.h"
#include "Core/Tensor.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Core/TensorBLAS.h"

class DeltaOperator {
public:
	DeltaOperator(const Tensorcd& Phi,
		const vector<SparseMatrixTreecd>& Hmats,
		const vector<SparseMatrixTreecd>& Holes,
		const SOPcd& sop, const Node& node)
		: holes_(Holes), node_(node) {

		npart_ = sop.size();
		vector<Tensorcd> Phis;
		for (size_t n = 0; n < sop.size(); ++n) {
			hPhi_.push_back(TreeFunctions::apply(Hmats[n], Phi, sop(n), node));
		}
	}

	~DeltaOperator() = default;

	size_t dim1() const {
		assert(!hPhi_.empty());
		return hPhi_[0].shape().lastBefore();
	}

	size_t dim2() const {
		return dim1();
	}

/*	Matrixcd operator*(const Matrixcd& r) const {

		cerr << "DeltaOperator::operator*: call product function instead.\n";
		getchar();
		Tensorcd R = toTensor(r);
		R.reshape(hPhi_[0].shape());

		Tensorcd Delta_R(R.shape());
		for (size_t l = 0; l < npart_; ++l) {
//			Matrixcd ml = hPhi_[l].dotProduct(R);
			Matrixcd ml = contractionBLAS(hPhi_[l], R, R.shape().lastIdx());
//			contractionBLAS(ml_, hPhi_[l], R, R.shape().lastIdx(), true);
			for (size_t o = 0; o < npart_; ++o) {
				const SparseMatrixTreecd& hole = holes_[npart_ * l + o];
				Matrixcd Mlo = hole[node_].transpose() * ml;
//				Delta_R += tensorMatrix(hPhi_[o], Mlo, node_.parentIdx());
				matrixTensorBLAS(Delta_R, Mlo.transpose(), hPhi_[o], node_.parentIdx(), false);
			}
		}
		Delta_R.reshape({r.dim1(), r.dim2()});
		return toMatrix(Delta_R);
	}
*/
	vector<Tensorcd> hPhi_;
	const vector<SparseMatrixTreecd>& holes_;
	size_t npart_;
	const Node& node_;
};

class DeltaMemory {
public:
	DeltaMemory(const Tensorcd& Phi,
		const vector<SparseMatrixTreecd>& Hmats,
		const SOPcd& sop, const Node& node) {

		Tensorcd hPhi = TreeFunctions::apply(Hmats[0], Phi, sop(0), node);
		size_t dim = hPhi.shape().lastDimension();
		ml_ = Matrixcd(dim, dim);
		workA = Tensorcd(hPhi.shape());
		workB = Tensorcd(hPhi.shape());
	}

	~DeltaMemory() = default;

	Matrixcd ml_;
	Tensorcd workA, workB;
};

Matrixcd product(const DeltaOperator& w, const Matrixcd& r, DeltaMemory *mem_ptr) {
	if (mem_ptr == nullptr) {
		cerr << "DeltaOperator error: mem is null.\n";
		exit(1);
	}
	DeltaMemory& mem = *mem_ptr;
	Tensorcd R(w.hPhi_[0].shape(), &r[0], false, false);

	Tensorcd Delta_R(R.shape());
	for (size_t l = 0; l < w.npart_; ++l) {
		contractionBLAS(mem.ml_, mem.workA, mem.workB, w.hPhi_[l], R, R.shape().lastIdx(), true);
		for (size_t o = 0; o < w.npart_; ++o) {
			const SparseMatrixTreecd& hole = w.holes_[w.npart_ * l + o];
			Matrixcd Mlo = hole[w.node_].transpose() * mem.ml_;
			matrixTensorBLAS(Delta_R, mem.workA, Mlo.transpose(), w.hPhi_[o], w.node_.parentIdx(), false);
		}
	}
	Delta_R.reshape({r.dim1(), r.dim2()});
	return moveToMatrix(Delta_R);
}


#endif //DELTAOPERATOR_H
