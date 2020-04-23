//
// Created by Hayley treir on 18/02/2020.
//
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeShape/TreeFactory.h"
#include "Util/RandomMatrices.h"

using namespace TreeFunctions;

void calculateDensity(MatrixTreecd& Rho_tree, const TensorTreecd& A_tree, const TensorTreecd& X_tree, const MatrixTreecd& S_tree,
	 const Matrixcd& S_label, const Tree& tree) {

	// top down iteration
	for (auto it = tree.rbegin(); it != tree.rend(); it++) {
		const Node& node = *it;

		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			ContractionLocal(Rho_tree, A_tree[parent], X_tree[parent], node, &S_tree);
		} else {
			// calculate delta S for top layer
			Rho_tree[node] = S_tree[node] - S_label;
		}
	}
}

void calculateGradient(TensorTreecd& Grad_tree, const TensorTreecd& A_tree, const TensorTreecd& X_tree, const MatrixTreecd& S_tree,
	const MatrixTreecd& Rho_tree, const Matrixcd& S_label, const Tree& tree) {

	for (auto it = tree.rbegin(); it != tree.rend(); it++) {
		const Node& node = *it;
		const TensorShape& tdim = node.shape();
		Tensorcd B = X_tree[node];
		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child_node = node.child(k);
				B = MatrixTensor(S_tree[child_node], B, k); // overwrite into A
			}
		}
		// grad = deltaS * B
		Grad_tree[node] = MatrixTensor(Rho_tree[node], B, tdim.lastIdx());
	}
}

void gradDescent() {

	/* Trees:
	 * 1. Weights tree, A_tree
	 * 2. Data tree, X_tree
	 * 3. Overlap tree, S_tree
	 * 4. Density tree, Rho_tree
	 * 5. Gradient tree, Grad_tree
	 */

	// generate tree
	Tree tree = TreeFactory::BalancedTree(3, 3, 2);

	// fill tree randomly
	mt19937 gen(1995);
	TensorTreecd A_tree(gen, tree);
	cout << "\n  => A_tree <= \n";
	A_tree.print(tree);

	// get data
	TensorTreecd X_tree(gen, tree, false);
	cout << "\n  => X <= \n";
	X_tree.print(tree);

	// calculate overlap, S (feedforward)
	using namespace TreeFunctions;
	MatrixTreecd S_tree = DotProduct(A_tree, X_tree, tree);
	cout << "\n  => S_tree <= \n";
	S_tree.print(tree);

	// generate randomised S_label matrix (will be passed in with correct labels)
	Matrixcd S_label = RandomMatrices::GUE(tree.TopNode().shape().lastDimension(), gen);
	cout << "\n  => S_label <= \n";
	S_label.print();

	// first pass down tree gets density matrices, Rho_tree
	MatrixTreecd Rho_tree(tree);
	calculateDensity(Rho_tree, A_tree, X_tree, S_tree, S_label, tree);
	cout << "\n  => Rho_tree <= \n";
	Rho_tree.print(tree);

	// calculate gradient
	TensorTreecd Grad_tree(tree);
	calculateGradient(Grad_tree, A_tree, X_tree, S_tree, Rho_tree, S_label, tree);
	cout << "\n  => Grad_tree <= \n";
	Grad_tree.print(tree);
}

int main() {
	gradDescent();
	return 0;
}