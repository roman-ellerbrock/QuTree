//
// Created by Grace Johnson on 1/30/20.
//
// TensorTree creation, population, basic operations, dot products, hole products
//

#include "TreeClasses/TensorTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeShape/TreeFactory.h"

// Demonstrate various ways to create a TensorTree object
TensorTreecd create_tensor_tree() {

    // 1. Create a TensorTree from a TensorTreeBasis (TTBasis) object
    size_t num_leaves = 4;
    size_t dim_leaves = 3;
    size_t dim_nodes = 2;
	Tree tree = TreeFactory::BalancedTree(num_leaves, dim_leaves, dim_nodes);
    tree.info();
    cout << "\nNo. total nodes: " << tree.nTotalNodes() << endl;
    cout << "No. logical nodes: " << tree.nNodes() << endl;
    cout << "No. leaves: " << tree.nLeaves() << endl;

    // Create tree and populate coefficients with zeroes
    TensorTreecd T(tree);
    cout << "\nT with zeroed coefficients :\n" << endl;
    T.print(tree); // print matrix elements for each node

    // Populate coefficients using a random number generator
    mt19937 gen(2468);
    TensorTreecd T_rand(tree);
	T_rand.FillRandom(gen, tree, false);
    cout << "\nT with random coefficients :\n" << endl;
    T_rand.print(tree);

    // 2. Read from an istream or file
    T_rand.Write("TT.dat");
    TensorTreecd T_rand2("TT.dat");

    return T_rand2;
}

void dot_product_tree() {
    cout << "\nTensor tree dot product:\n" << endl;

    // A dot product between two tensor trees results in a factor matrix tree
	Tree tree = TreeFactory::BalancedTree(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(gen, tree, false);
    TensorTreecd Chi(Psi);

    /// Construct, allocate and calculate a factor matrix tree
    MatrixTreecd w = MatrixTreeFunctions::DotProduct(Psi, Chi, tree);
    w.print();

    return;
}

void hole_product_tree() {
    cout << "\nTensor tree hole product:\n" << endl;

    // A hole product between two tensor trees results in a hole matrix tree
	Tree tree = TreeFactory::BalancedTree(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(gen, tree, false);
    TensorTreecd Chi(Psi);
//    FactorMatrixTreecd w(Psi, Chi, tree);
    MatrixTreecd w = MatrixTreeFunctions::DotProduct(Psi, Chi, tree);

    /// Construct, allocate and calculate a hole matrix tree
//    HoleMatrixTreecd rho(Psi, Chi, w, tree);
    MatrixTreecd rho = MatrixTreeFunctions::Contraction(Psi, Chi,  w, tree);
    rho.print();

}

#include "TreeShape/TreeFactory.h"

void tree_examples() {
	Tree tree = TreeFactory::BalancedTree(3, 3, 2);
	TensorTreecd Psi(tree);
	Psi.print(tree);

	mt19937 gen(1995);
	Psi.FillRandom(gen, tree);
	Psi.print(tree);

	cout << "A MatrixTree from an overlap:\n";
	using namespace MatrixTreeFunctions;
	MatrixTreecd W = DotProduct( Psi, Psi, tree);
	W.print(tree);

	cout << "Density Matrix:\n";
	MatrixTreecd Rho = Contraction(Psi, tree, true);
	Rho.print(tree);

	cout << "<Psi|Chi> MatrixTree:\n";
	TensorTreecd Chi(gen, tree, false);
	MatrixTreecd S = DotProduct(Psi, Chi, tree);
	S.print(tree);

	cout << "Contractions of Psi and Chi:\n";
	MatrixTree C = Contraction(Psi, Chi, S, tree);
	C.print(tree);
}

int main() {
	tree_examples();

/*    TensorTreecd T = create_tensor_tree();
    TensorTreecd T2(T);
    dot_product_tree();
    hole_product_tree();

    Tree tree(7, 2, 2);
	Matrixcd X(2, 2);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, 3);

	mt19937 gen(2020);
	TensorTreecd Psi(tree, gen);

	SparseMatrixTreecd hmats(M, tree);
	SparseMatrixTreeFunctions::Represent(hmats, M, Psi, Psi, tree);
	hmats.print();
	// TODO: sparse versions of dot product and hole product (initialized with operators)
*/

}