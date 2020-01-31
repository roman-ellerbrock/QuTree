//
// Created by Grace Johnson on 1/30/20.
//
// TensorTree creation, population, basic operations, dot products, hole products
//

#include "Tree/TensorTree.h"

// Demonstrate various ways to create a TensorTree object
TensorTreecd create_tensor_tree() {

    // 1. Create a TensorTree from a TensorTreeBasis (TTBasis) object
    size_t num_leaves = 4;
    size_t dim_leaves = 3;
    size_t dim_nodes = 2;
    TTBasis basis(num_leaves, dim_leaves, dim_nodes); // Creates a balanced tree by default
    basis.info();
    cout << "\nNo. total nodes: " << basis.nTotalNodes() << endl;
    cout << "No. logical nodes: " << basis.nNodes() << endl;
    cout << "No. leaves: " << basis.nLeaves() << endl;

    // Create tree and populate coefficients with zeroes
    TensorTreecd T(basis);
    cout << "\nT with zeroed coefficients :\n" << endl;
    T.print(basis); // print matrix elements for each node

    // Populate coefficients using a random number generator
    mt19937 gen(2468);
    TensorTreecd T_rand(basis);
    T_rand.Generate(basis, gen, false);
    cout << "\nT with random coefficients :\n" << endl;
    T_rand.print(basis);

    // 2. Read from an istream or file
    T_rand.Write("TT.dat");
    TensorTreecd T_rand2("TT.dat");

    return T_rand2;
}

void dot_product_tree(TensorTreecd A, TensorTreecd B) {
    // TODO:
    return;
}

void hole_product_tree(TensorTreecd A, TensorTreecd B) {
    // TODO:
    return;
}

int main() {

    TensorTreecd T = create_tensor_tree();
    TensorTreecd T2(T);
    dot_product_tree(T, T2);
    hole_product_tree(T, T2);

}