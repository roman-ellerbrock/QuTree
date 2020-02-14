//
// Created by Grace Johnson on 1/30/20.
//
// TensorTree creation, population, basic operations, dot products, hole products
//

#include "Tree/TensorTree.h"
#include "Tree/MatrixTreeFunctions.h"
#include "Tree/SparseMatrixTreeFunctions.h"

// Demonstrate various ways to create a TensorTree object
TensorTreecd create_tensor_tree() {

    // 1. Create a TensorTree from a TensorTreeBasis (TTBasis) object
    size_t num_leaves = 4;
    size_t dim_leaves = 3;
    size_t dim_nodes = 2;
    TTBasis tree(num_leaves, dim_leaves, dim_nodes); // Creates a balanced tree by default
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
    T_rand.Generate(tree, gen, false);
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
    TTBasis tree(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(tree, gen, false);
    TensorTreecd Chi(Psi);

    /// Construct, allocate and calculate a factor matrix tree
    MatrixTreecd w = MatrixTreeFunctions::DotProduct(Psi, Chi, tree);
    w.print();

    return;
}

void hole_product_tree() {
    cout << "\nTensor tree hole product:\n" << endl;

    // A hole product between two tensor trees results in a hole matrix tree
    TTBasis tree(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(tree, gen, false);
    TensorTreecd Chi(Psi);
//    FactorMatrixTreecd w(Psi, Chi, tree);
    MatrixTreecd w = MatrixTreeFunctions::DotProduct(Psi, Chi, tree);

    /// Construct, allocate and calculate a hole matrix tree
//    HoleMatrixTreecd rho(Psi, Chi, w, tree);
    MatrixTreecd rho = MatrixTreeFunctions::Contraction(Psi, Chi,  w, tree);
    rho.print();

}

int main() {

    TensorTreecd T = create_tensor_tree();
    TensorTreecd T2(T);
    dot_product_tree();
    hole_product_tree();

    TTBasis tree(7, 2, 2);
	FactorMatrixcd X(2, 1);
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


}