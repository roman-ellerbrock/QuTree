//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "DVR/cdvr_functions.h"
#include "Core/Tensor_Extension.h"
#include "Core/TensorBLAS.h"

namespace cdvr_functions {

	/**
	 * Notes on notation:
	 * Cdown: Top-down normalized, reads C^{p(p\circ k)} in equations
	 * Cup: Bottom-up normalized, reads C^{p(p\circ 0)} in equations
	 */

	/**
	 * Fill a vector with grid points corresponding to intex "idx".
	 * @param X
	 * @param idx
	 * @param grids
	 * @param node
	 */

	void fillX(Vectord& X, size_t idx, const TreeGrids& grids, const Node& node) {
		/// Fill full-dimensional grid points from all neighboring node-grids
		for (size_t k = 0; k < X.dim(); ++k) {
			const SparseVectorTreed& grid = grids[k];
			if (grid.isActive(node)) {
				const Vectord& localgrid = grid[node];
				X(k) = localgrid(idx);
			}
		}
	}

	void fillXNode(Vectord& X, vector<size_t> idx, const TreeGrids& grids,
		const TreeGrids& holegrids, const Node& node) {
		/// \brief Fill X for a node-grid
		assert(grids.size() == holegrids.size());
		assert(X.dim() == grids.size());
		if (node.isBottomlayer()) {
			const Leaf& leaf = node.getLeaf();
			const auto& g = leaf.interface();
			assert(g.hasDVR());
			const Vectord& x = g.getX();
			X(leaf.mode()) = x(idx.front());
		} else {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				fillX(X, idx[k], grids, child);
			}
		}

		fillX(X, idx.back(), holegrids, node);
	}

	void fillXEdge(Vectord& X, vector<size_t> idx, const TreeGrids& grids,
		const TreeGrids& holegrids, const Node& node) {
		/// \brief Fill X for an edge-grid
		assert(grids.size() == holegrids.size());
		assert(X.dim() == grids.size());
		assert(idx.size() == 2);
		fillX(X, idx.front(), grids, node);
		fillX(X, idx.back(), holegrids, node);
	}

	Tensorcd buildE(const Tensorcd& subDeltaV,
		const Tensorcd& D, const Node& node, size_t k) {

		const TensorShape& shape = node.shape();
		size_t dimc = shape[k];
		size_t dimn = shape[node.nChildren()];
		TensorShape eshape({dimc, dimc, dimn, dimn});
		Tensorcd E(eshape);

		size_t dimc2 = dimc * dimc;
		size_t dimc3 = dimc * dimc * dimc;
		// eshape = {dimc, dimc, dimn, dimn}
		// subdeltaVshape = {dimc, dimc, dimc, dimc}
		// dshape = {dimc, dimn, dimc, dimn}
		for (size_t i0 = 0; i0 < dimc; ++i0) {
			for (size_t i1 = 0; i1 < dimc; ++i1) {
				for (size_t i2 = 0; i2 < dimn; ++i2) { // dimn
					for (size_t i3 = 0; i3 < dimn; ++i3) { // dimn
						for (size_t l1 = 0; l1 < dimc; ++l1) {
							for (size_t l2 = 0; l2 < dimc; ++l2) {
								// vidx = {idx[0], idx[1], l1, l2};
								const size_t vidx = i0 + i1 * dimc + l1 * dimc2 + l2 * dimc3;
								// didx = {l2, idx[3], l1, idx[2]};
								const size_t didx = l2 + i3 * dimc + l1 * dimc * dimn + i2 * dimc2 * dimn;
								// indexMapping(idx, J, eshape);
								const size_t J = i0 + i1 * dimc + i2 * dimc2 + i3 * dimc2 * dimn;
								E(J) += subDeltaV(vidx) * D(didx);
							}
						}
					}
				}
			}
		}

		return E;
	}

	void contractEF(Tensorcd& deltaV, const Tensorcd& E, const Tensorcd& F,
		const Node& node, size_t k) {

		const Node& child = node.child(k);
		const TensorShape& shape = node.shape();
		size_t dimc = shape[k];
		const TensorShape& vshape = deltaV.shape();
		size_t dimn = vshape[0];
		size_t dimc2 = dimc * dimc;
		size_t dimn2 = dimn * dimn;
		size_t dimn3 = dimn * dimn * dimn;
		for (size_t i0 = 0; i0 < dimn; ++i0) {
			for (size_t i1 = 0; i1 < dimn; ++i1) {
				for (size_t i2 = 0; i2 < dimn; ++i2) {
					for (size_t i3 = 0; i3 < dimn; ++i3) {
						for (size_t l1 = 0; l1 < dimc; ++l1) {
							for (size_t l2 = 0; l2 < dimc; ++l2) {
//								size_t I = indexMapping({i0, i1, i2, i3}, vshape);
								size_t I = i0 + i1 * dimn + i2 * dimn2 + i3 * dimn3;
								size_t fidx = l1 + i0 * dimc + l2 * dimc * dimn + i1 * dimc2 * dimn;
								size_t eidx = l1 + l2 * dimc + i2 * dimc2 + i3 * dimc2 * dimn;
								deltaV(I) += F(fidx) * E(eidx);
//								deltaV(I) += F({l1, i0, l2, i1}) * E({l1, l2, i2, i3});
							}
						}
					}
				}
			}
		}
	}

	void deltaEdgeCorrection(Tensorcd& deltaV, const Tensorcd& Cup,
		const TensorTreecd& Cdown, const DeltaVTree& DeltaVs,
		const Node& node) {

		assert(!node.isBottomlayer());
		const TensorShape& shape = node.shape();
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			const Tensorcd& SubDeltaV = DeltaVs[child];

			auto D = Tensor_Extension::doubleHoleContraction(Cdown[child], Cup, k, node.nChildren());
			auto F = Tensor_Extension::doubleHoleContraction(Cup, Cdown[child], k, node.nChildren());
			Tensorcd E = buildE(SubDeltaV, D, node, k);
			contractEF(deltaV, E, F, node, k);
		}
	}

	void calculateDeltaEdgeLocal(Tensorcd& deltaV, const Tensorcd& Cup,
		const Tensorcd& Vnode, const Matrixd& Vedge, const Node& node) {
		/**
		 * \brief Calculate deltaV-tensor for bottomlayer node
		 * @param deltaV "matrix"-representation of operator
		 * @param C Wavefunction coefficients at node
		 * @param Vnode Potential evaluated at node-DVR
		 * @param Vedge Potential evaluated at edge-DVR
		 * @param node Node at which the edge is pointing to (std convention for edges in trees)
		 *
		 * There are no underlying deltaV matrices at node. Therefore, a simplified expression.
		 */

		const TensorShape& shape = Vnode.shape();
		const TensorShape& Cshape = Cup.shape();
		size_t dim = shape.lastDimension();
		size_t dimbef = shape.lastBefore();
		assert(shape.totalDimension() == Cshape.totalDimension());

		/// Sanity checks
		assert(!node.isToplayer());
		assert(deltaV.shape().totalDimension() == pow(shape.lastDimension(), 4));

		deltaV.zero();
		vector<size_t> idxs(4);
		for (size_t l1 = 0; l1 < dim; ++l1) {
			for (size_t i0 = 0; i0 < dim; ++i0) {
				for (size_t m1 = 0; m1 < dim; ++m1) {
					idxs = {l1, i0, m1, i0};
					size_t J = indexMapping(idxs, deltaV.shape());
					for (size_t Ibef = 0; Ibef < dimbef; Ibef++) {
						deltaV(J) += conj(Cup(Ibef, l1)) * Vnode(Ibef, i0) * Cup(Ibef, m1);
					}
				}
			}
		}

		/// Substract Vedge from diagonal
		for (size_t l0 = 0; l0 < dim; ++l0) {
			for (size_t l1 = 0; l1 < dim; ++l1) {
				idxs = {l0, l1, l0, l1};
				size_t L = indexMapping(idxs, deltaV.shape());
				deltaV(L) -= Vedge(l0, l1);
			}
		}
	}

	void calculateDeltaVs(DeltaVTree& deltaVs,
		const TensorTreecd& Cup, const TensorTreecd& Cdown,
		const TensorTreecd& Vnodes, const MatrixTreed& Vedges,
		const Tree& tree) {

		for (const Node& node: tree) {
			if (node.isBottomlayer()) {
				calculateDeltaEdgeLocal(deltaVs[node], Cup[node],
					Vnodes[node], Vedges[node], node);
			} else if (!node.isToplayer()) {
				calculateDeltaEdgeLocal(deltaVs[node], Cup[node],
					Vnodes[node], Vedges[node], node);
				deltaEdgeCorrection(deltaVs[node], Cup[node], Cdown,
					deltaVs, node);
			}
		}
	}

	Matrixcd applyDeltaV(const Tensorcd& deltaV, const Matrixcd& x,
		size_t dim) {

		const TensorShape& shape = deltaV.shape();
		Matrixcd y(dim, dim);

		size_t dim2 = dim * dim;
		size_t dim3 = dim * dim * dim;
		for (size_t l0 = 0; l0 < dim; ++l0) {
			for (size_t l1 = 0; l1 < dim; ++l1) {
				for (size_t l2 = 0; l2 < dim; ++l2) {
					for (size_t l3 = 0; l3 < dim; ++l3) {
						size_t L = l0 + l1 * dim + l2 * dim2 + l3 * dim3;
						y(l0, l1) += deltaV(L) * x(l3, l2);
					}
				}
			}
		}

		return y;
	}

	void applyCorrection(Tensorcd& VPhi, const Tensorcd& Phi, const Tensorcd& C,
		const Tensorcd& deltaV, const Node& child, const WorkMemorycd& mem) {

		/// I: Contract over
		size_t k = child.childIdx();
		const TensorShape& shape = deltaV.shape();
		assert(C.shape().totalDimension() == Phi.shape().totalDimension());

		Matrixcd x = contractionBLAS(C, Phi, k);

		/// II: apply DeltaV
		Matrixcd y = applyDeltaV(deltaV, x, Phi.shape()[k]);

		/// III: M * C
		VPhi += matrixTensorBLAS(y, C, k);
	}

	void apply(Tensorcd& VXi, const Tensorcd& Xi, const Tensorcd& V,
		const TensorTreecd& Cdown, const DeltaVTree& deltaVs,
		const Node& node, const WorkMemorycd& mem) {

		/// TODO: V diagonal in last idx for toplayer
		VXi = productElementwise(V, Xi);

		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				/// VXi += deltaV_k
				applyCorrection(VXi, Xi, Cdown[child], deltaVs[child], child, mem);
			}
		}
	}

/*	TensorTreecd Apply(const Wavefunction& Psi, const TensorTreecd& Cdown,
		const TensorTreecd& V, const DeltaVTree& DeltaVs, const Tree& tree) {

		Wavefunction VPsi = Psi;
		for (const Node& node : tree) {
			VPsi[node] = Apply(Psi[node], V[node], Cdown, DeltaVs, node);
		}

		return VPsi;
	}*/
}

