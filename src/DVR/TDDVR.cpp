//
// Created by Roman Ellerbrock on 3/9/20.
//

#include <Core/TensorBLAS.h>
#include "DVR/TDDVR.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Util/WeightedSimultaneousDiagonalization.h"
#include "TreeClasses/SymMatrixTreeFunctions.h"

vector<Matrixcd> getXs(const vector<SparseMatrixTreecd>& Xs, const Node& node) {
	vector<Matrixcd> xs;
	for (const auto& xtree: Xs) {
		if (xtree.isActive(node)) {
			xs.push_back(xtree[node]);
		}
	}
	return xs;
}

void setGrids(TreeGrids& gridtrees, vector<Vectord> grids, const Node& node) {
	size_t k = 0;
	for (auto& grid: gridtrees) {
		if (grid.isActive(node)) {
			grid[node] = grids[k];
			k++;
		}
	}
}

vector<double> calculateShift(const vector<Matrixcd>& xs, const Matrixcd& w) {
	vector<double> shifts;
	for (const auto& x: xs) {
		double s = abs(x.trace());
		shifts.push_back(s);
	}
	return shifts;
}

void shift(vector<Matrixcd>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Matrixcd& x = xs[i];
		for (size_t j = 0; j < x.dim1(); ++j) {
			x(j, j) -= s;
		}
	}
}

void shiftBack(vector<Vectord>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Vectord& x = xs[i];
		for (size_t j = 0; j < x.dim(); ++j) {
			x(j) += s;
		}
	}
}

void LayerGrid(TreeGrids& grids, Matrixcd& trafo,
	const vector<SparseMatrixTreecd>& Xs,
	Matrixcd w, const Node& node) {
	assert(grids.size() == Xs.size());
	auto xs = getXs(Xs, node);
	assert(!xs.empty());

	if (xs.size() == 1) {
		SpectralDecompositioncd diags = diagonalize(xs.front());
		trafo = diags.first;
		trafo = trafo.adjoint();
		setGrids(grids, {diags.second}, node);
	} else {
		w = regularize(w, 1e-9);

		auto shifts = calculateShift(xs, w);
		shift(xs, shifts);

		auto diags = WeightedSimultaneousDiagonalization::calculate(xs, w, 1e-5);

		shiftBack(diags.second, shifts);

		setGrids(grids, diags.second, node);
		trafo = diags.first;
		trafo = trafo.adjoint();
	}
}

void updateGrids(TreeGrids& grids, MatrixTreecd& trafo, const vector<SparseMatrixTreecd>& Xs,
	const MatrixTreecd& rho, const Tree& tree) {
	for (const Node& node: tree) {
		if (node.isToplayer()) { continue; }
		LayerGrid(grids, trafo[node], Xs, rho[node], node);
	}
}

void TDDVR::update(const Wavefunction& Psi, const Tree& tree) {
	/// Calculate density matrix
	TreeFunctions::contraction(rho_, Psi, tree, true);

	/// Calculate X-Matrices
	Xs_.Update(Psi, tree);

	/// Build standard grid
	updateGrids(grids_, trafo_, Xs_.mats_, rho_, tree);

	/// Build hole grid
	updateGrids(hole_grids_, hole_trafo_, Xs_.holes_, rho_, tree);
}

void layerGrid(TreeGrids& grids, Matrixcd& trafo,
	vector<Matrixcd> xs,
	Matrixcd w, const Node& node) {

	assert(grids.size() == Xs.size());
	assert(!xs.empty());

	if (xs.size() == 1) {
		SpectralDecompositioncd diags = diagonalize(xs.front());
		trafo = diags.first;
		trafo = trafo.adjoint();
		setGrids(grids, {diags.second}, node);
	} else {
		w = regularize(w, 1e-9);

		auto shifts = calculateShift(xs, w);
		shift(xs, shifts);

		auto diags = WeightedSimultaneousDiagonalization::calculate(xs, w, 1e-5);

		shiftBack(diags.second, shifts);

		setGrids(grids, diags.second, node);
		trafo = diags.first;
		trafo = trafo.adjoint();
	}
}

vector<Matrixcd> getX(const SymXMatrixTrees& Xs,
	const Node& node, bool up) {
	vector<Matrixcd> xs;
	if (up) {
		for (const auto& xsym: Xs.xmat_) {
			const auto& xtree = xsym.up();
			if (xtree.isActive(node)) {
				xs.push_back(xtree[node]);
			}
		}
	} else {
		for (const auto& xsym: Xs.xmat_) {
			const auto& xtree = xsym.down();
			if (xtree.isActive(node)) {
				xs.push_back(xtree[node]);
			}
		}
	}
	return xs;
}

void updateGrids(SymTreeGrid& sgrids, SymMatrixTree& u, const SymXMatrixTrees& xsym,
	const SymMatrixTree& rho, const Tree& tree) {

	for (const Node& node: tree) {
		auto xs = getX(xsym, node, true);
		layerGrid(sgrids.up(), u.up()[node], xs, rho.down()[node], node);
	}

	for (const Node& node: tree) {
		auto xs = getX(xsym, node, false);
		layerGrid(sgrids.down(), u.down()[node], xs, rho.up()[node], node);
	}
}

void TDDVR::update(const SymTensorTree& Psi, const Tree& tree) {

	/// Calculate sym density-Matrices
	auto srho = TreeFunctions::weightedContraction(Psi, Psi, tree);

	/// Calculate X-Matrices
	symx_.update(Psi, tree);

	/// Evaluate Grids
	updateGrids(sgrids_, strafo_, symx_, srho, tree);
}

void TDDVR::NodeTransformation(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// @TODO: make in/out independent
	/// Transform underlying A-coefficient
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (!inverse) {
				Phi = matrixTensorBLAS(trafo_[child], Phi, k);
			} else {
				Phi = matrixTensorBLAS(trafo_[child].adjoint(), Phi, k);
			}
		}
	}

	/// Transform state
	if (!node.isToplayer()) {
		if (!inverse) {
			Phi = matrixTensorBLAS(hole_trafo_[node], Phi, node.nChildren());
		} else {
			Phi = matrixTensorBLAS(hole_trafo_[node].adjoint(), Phi, node.nChildren());
		}
	}
}

void TDDVR::downTransformation(Tensorcd& Phi, const Node& node, bool inverse) const {
	assert(!node.isToplayer());
	const Node& parent = node.parent();
	for (size_t k = 0; k < parent.nChildren(); ++k) {
		const Node& child = parent.child(k);
		Matrixcd U = inverse ? strafo_.down()[child].adjoint() : strafo_.down()[child];
		Phi = matrixTensor(U, Phi, k);
	}

	if (!parent.isToplayer()) {
		Matrixcd U = inverse ? strafo_.down()[node].adjoint() : strafo_.down()[node];
		Phi = matrixTensor(U, Phi, node.parentIdx());
	}
}

void TDDVR::EdgeTransformation(Matrixcd& B_inv, const Edge& edge, bool inverse) const {
	if (!inverse) {
		B_inv = trafo_[edge].transpose().adjoint() * B_inv;
		B_inv = B_inv * hole_trafo_[edge].adjoint();
	} else {
		B_inv = trafo_[edge].transpose() * B_inv;
		B_inv = B_inv * hole_trafo_[edge];
	}
}

void TDDVR::upTransformation(Tensorcd& Phi, const Node& node, bool inverse) const {
	assert(!node.isToplayer());
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
//			Matrixcd u = inverse ? strafo_.up()[child].adjoint() : strafo_.up()[child];
			Matrixcd u = inverse ? trafo_[child].adjoint() : trafo_[child];
			Phi = matrixTensor(u, Phi, k);
		}
	}

	if (!node.isToplayer()) {
		Matrixcd U = inverse ? hole_trafo_[node].adjoint() : hole_trafo_[node];
//		Matrixcd U = inverse ? strafo_.down()[node].adjoint() : strafo_.down()[node];
		Phi = matrixTensor(U, Phi, node.parentIdx());
	}
}

void TDDVR::downTransformation(SymTensorTree& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node: tree) {
		if (node.isToplayer()) { continue; }
		downTransformation(Psi.down_[node], node, inverse);
	}
}

void TDDVR::upTransformation(SymTensorTree& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node: tree) {
		upTransformation(Psi.up_[node], node, inverse);
	}
}

void TDDVR::NodeTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node: tree) {
		NodeTransformation(Psi[node], node, inverse);
	}
}

void TDDVR::EdgeTransformation(MatrixTreecd& B_inv, const Tree& tree, bool inverse) const {
	for (const Edge& edge: tree.edges()) {
		EdgeTransformation(B_inv[edge], edge, inverse);
	}
}

void TDDVR::GridTransformation(MatrixTensorTree& Psi, const Tree& tree, bool inverse) const {
	NodeTransformation(Psi.nodes(), tree, inverse);
	EdgeTransformation(Psi.edges(), tree, inverse);
}

void TDDVR::GridTransformation(SymTensorTree& Psi, const Tree& tree, bool inverse) const {
	upTransformation(Psi, tree, inverse);
	downTransformation(Psi, tree, inverse);
}

void TDDVR::print(const Tree& tree) const {
	cout << "TDDVR: " << endl;
	cout << "Grids:" << endl;
	for (const Node& node: tree) {
		if (!node.isToplayer()) {
			size_t dim = trafo_[node].dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid: grids_) {
					if (grid.isActive(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}
	cout << "Hole grids:" << endl;
	for (const Node& node: tree) {
		if (!node.isToplayer()) {
			size_t dim = trafo_[node].dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid: hole_grids_) {
					if (grid.isActive(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}
}

/*void TDDVR::GridTransformationLocal(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// Transform underlying A-coefficient
	for (size_t k = 0; k < node.nChildren(); ++k) {
		const Node& child = node.child(k);
		if (!inverse) {
			Phi = matrixTensorBLAS(trafo_[child], Phi, k);
		} else {
			Phi = matrixTensorBLAS(trafo_[child].adjoint(), Phi, k);
		}
	}

	/// Transform state
	if (!inverse) {
//		Phi = tensorMatrix(Phi, hole_trafo_[node], node.parentIdx());
		Phi = matrixTensorBLAS(hole_trafo_[node].transpose(), Phi, node.parentIdx());
	} else {
//		Phi = tensorMatrix(Phi, hole_trafo_[node].adjoint(), node.parentIdx());
		Phi = matrixTensorBLAS(hole_trafo_[node].adjoint().transpose(), Phi, node.parentIdx());
	}
}*/

/*void TDDVR::GridTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node : tree) {
		GridTransformationLocal(Psi[node], node, inverse);
	}
}*/



