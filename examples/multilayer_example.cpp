//
// Created by Roman Ellerbrock on 1/22/21.
//

#include "TreeClasses/TensorTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/SpectralDecompositionTree.h"

namespace examples {

	/**
	 * Rationale:
	 * (a) This is an example that guides though Ref. 1 and
	 * 	   shows how the Equations are implemented in QuTree.
	 * (b)
	 *
	 * References:
	 * [1] U. Manthe, J. Chem. Phys. 128, 164116 (2008)
	 */

	/// Small helper to generate X- and 1-Gate
	LeafMatrixcd xGate() {
		Matrixcd X(2, 2);
		X(0, 1) = 1.;
		X(1, 0) = 1.;
		return LeafMatrixcd(X);
	}
	LeafMatrixcd oneGate() {
		Matrixcd I(2, 2);
		I(0, 0) = 1.;
		I(1, 1) = 1.;
		return LeafMatrixcd(I);
	}

	void Manthe2008Paper() {
		/// We start by generating the tree shown in Fig. 1 f).
		/// We choose f=12 leaves, N=2, n=3 and a=1 ("a" is can not be modified in this function).
		Tree tree = TreeFactory::BalancedTree(12, 2, 3);


		/// === Wave Function Generation & Basic Operations (mostly Sec. II A.) ===
		/// Generate a random wave function as introduced in Sec. II A.
		mt19937 gen(92349);
		bool HartreeProduct = false; /// true: generate a Hartree Product, false: generate random tensors
		TensorTreecd Psi(gen, tree, HartreeProduct);
		/// For demonstration purposes, we generate a second wavefunction
		TensorTreecd Chi(gen, tree, HartreeProduct);

		/// This calculates the overlap-matrices resulting from the overlap of
		/// two wave functions <Psi|Chi>. Note: this is not displayed in Ref. 1.
		MatrixTreecd S = TreeFunctions::DotProduct(Psi, Chi, tree);
		/// The overlap of the two wave functions is the overlap matrix at the toplayer
		Matrixcd braket_mat = S[tree.TopNode()];
		/// Note:
		/// (a) <Psi|Chi> is in general a matrix, if the state-averaged formulation is used.
		///     If, like here, we calculate the overlap of single wave functions, it is a 1x1-Matrix.
		/// (b) This is not shown in Ref. 1.
		/// (c) Psi (and Chi) are generated with orthonormal SPF basis sets, therefore S=Dot(Psi, Psi, tree)
		///     would only consist of identity matrices.

		/// Reduced density matrices as displayed in Eq. (19) of Ref. 1 can be calculated by:
		MatrixTreecd rho = TreeFunctions::Contraction(Psi, tree, true); /// true: SPF basis is orthogonal

		/// The wave function is transformed to natural orbitals via
		CanonicalTransformation(Psi, tree);
		/// ...after re-calculating density matrices, they are diagonal.
		rho = TreeFunctions::Contraction(Psi, tree, true);


		/// === Access information in Tree-Objects ===
		/// Tree-structured objects are maps from Nodes to node-attributes (or the same for edges).
		/// To access information at a node, one typically swipes through the tree via
		for (const Node& node : tree) {
			Tensorcd& Phi = Psi[node];
			Matrixcd& s = S[node];
			/// and so on...
			/// Adjacent nodes are acced via
			if (!node.isToplayer()) {
				const Node& parent = node.parent(); /// get the parent
				Tensorcd& Phi_parent = Psi[parent];
			}
			if (!node.isBottomlayer()) {
				const Node& child = node.child(0); /// get the first child
				Tensorcd& Phi_child = Psi[child];
			}
		}
		/// The standard iterator performs a bottom-up swipe; there is also a reverse iterator
		/// to perform a top-down swipe.


		/// === Operator Representation (Sec. II C.) ===
		/// Here we show how to create an operator and represent it.
		/// We restrict the demonstration to a single summand of a SOP operator.
		/// A simple operator acting on one leaf is called a LeafOperator.
		/// A product of LeafOperators is called a MultiLeafOperators (MLO in short).
		/// A sum of MLOs is a sum-of-product operator.
		MLOcd product_operator(xGate(), 0);
		product_operator.push_back(xGate(), 3);/// X acting on qubit 0 and 3, X(0)*X(3)
		/// This is a single summand in Eq. (28) of Ref. 1.

		/// The matrix elements of this operator (Eqs. (32) & (33)) are obtained by
		SparseMatrixTreecd hmat = TreeFunctions::Represent(product_operator, Psi, Psi, tree);
		/// h-matrices are only required for a subset of nodes and therefore stored as a SparseMatrixTree.
		/// The contraction of h-matrices (hole matrices, Eqs. (36)) are obtained by
		SparseMatrixTreecd hhole = TreeFunctions::Contraction(Psi, hmat, tree);

		/// There is a lot more about operators, including overloadings, functions of operators, etc.
		/// However, I will not go into more detail here.
	}

}

int main() {
	examples::Manthe2008Paper();
	return 0;
}
