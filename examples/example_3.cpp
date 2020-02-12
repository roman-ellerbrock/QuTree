//
// Created by Grace Johnson on 1/30/20.
//
// TensorTree creation, population, basic operations, dot products, hole products
//

#include "Tree/TensorTree.h"
#include "Tree/FactorMatrixTree.h"
#include "Tree/HoleMatrixTree.h"

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

void dot_product_tree() {
    cout << "\nTensor tree dot product:\n" << endl;

    // A dot product between two tensor trees results in a factor matrix tree
    TTBasis basis(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(basis, gen, false);
    TensorTreecd Chi(Psi);

    /// Construct, allocate and calculate a factor matrix tree
    FactorMatrixTreecd w(Psi, Chi, basis);
    w.print();

    return;
}

void hole_product_tree() {
    cout << "\nTensor tree hole product:\n" << endl;

    // A hole product between two tensor trees results in a hole matrix tree
    TTBasis basis(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(basis, gen, false);
    TensorTreecd Chi(Psi);
    FactorMatrixTreecd w(Psi, Chi, basis);

    /// Construct, allocate and calculate a hole matrix tree
    HoleMatrixTreecd rho(Psi, Chi, w, basis);
    rho.print();

    return;
}

#include "MatrixTree.h"
#include "MatrixTreeFunctions.h"

int main() {

    TensorTreecd T = create_tensor_tree();
    TensorTreecd T2(T);
    dot_product_tree();
    hole_product_tree();

    TTBasis basis(4, 3, 2);
    mt19937 gen(2020);
    TensorTreecd Psi(basis, gen, false);
	TensorTreecd Chi(basis, gen, false);
	FactorMatrixTreecd S(Psi, Chi, basis);

    MatrixTreecd S2(basis);
    MatrixTreeFunctions::DotProduct(S2, Psi, Chi, basis);
	S.print(basis);
	S2.print(basis);

	HoleMatrixTreecd Rho(Psi, Chi, S, basis);

	MatrixTreecd Rho2(basis);
	MatrixTreeFunctions::Contraction(Rho2, Psi, Chi, S2, basis);

	cout << "Rho:\n";
	Rho.print(basis);
	cout << "Rho2:\n";
	Rho2.print(basis);

	// TODO: sparse versions of dot product and hole product (initialized with operators)

}