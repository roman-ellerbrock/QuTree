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
    MatrixTreecd w = TreeFunctions::dotProduct(Psi, Chi, tree);
    w.print();
}

void hole_product_tree() {
    cout << "\nTensor tree hole product:\n" << endl;

    // A hole product between two tensor trees results in a hole matrix tree
	Tree tree = TreeFactory::BalancedTree(4, 3, 2);
    mt19937 gen(2468);
    TensorTreecd Psi(gen, tree, false);
    TensorTreecd Chi(Psi);
//    FactorMatrixTreecd w(Psi, Chi, tree);
    MatrixTreecd w = TreeFunctions::dotProduct(Psi, Chi, tree);

    /// Construct, allocate and calculate a hole matrix tree
//    HoleMatrixTreecd rho(Psi, Chi, w, tree);
    MatrixTreecd rho = TreeFunctions::contraction(Psi, Chi, w, tree);
    rho.print();

}

void tree_examples() {
	Tree tree = TreeFactory::BalancedTree(12, 5, 2);
	TensorTreecd Psi(tree);
	Psi.print(tree);

	mt19937 gen(1995);
	Psi.FillRandom(gen, tree);
	Psi.print(tree);

	cout << "A MatrixTree from an overlap:\n";
	using namespace TreeFunctions;
	MatrixTreecd W = dotProduct(Psi, Psi, tree);
	W.print(tree);

	cout << "Density Matrix:\n";
	MatrixTreecd Rho = contraction(Psi, tree, true);
	Rho.print(tree);

	cout << "<Psi|Chi> MatrixTree:\n";
	TensorTreecd Chi(gen, tree, false);
	MatrixTreecd S = dotProduct(Psi, Chi, tree);
	S.print(tree);

	cout << "Contractions of Psi and Chi:\n";
	MatrixTreecd C = contraction(Psi, Chi, S, tree);
	C.print(tree);

	cout << "Add two tensor trees:\n";
	auto Eta = Psi + Chi;
	Eta.print(tree);

	cout << "Substract two tensor trees:\n";
	auto Beta = Psi - Chi;
	Beta.print(tree);

	cout << "Substract two tensor trees:\n";
	complex<double> coeff = 2.;
	auto Gamma = coeff * Beta;
	Beta.print(tree);
}

void paper_1() {
	/// generate Tree and TensorTree
	Tree tree = TreeFactory::BalancedTree(12, 5, 2);
	TensorTreecd Psi(tree);

	/// Operations on Tensor Trees
	using namespace TreeFunctions;
	MatrixTreecd W = dotProduct(Psi, Psi, tree); /// <Psi|Psi>_p
	MatrixTreecd C = contraction(Psi, Psi, W, tree); /// <Psi|Psi>_(p)
	MatrixTreecd Rho = contraction(Psi, tree, true); /// Assume orthogonal basis

	/// Represent Operators
	auto x = &LeafInterface::applyX;
	MultiLeafOperator<complex<double>> M(x, 0);
	auto Mrep = TreeFunctions::represent(M, Psi, tree); /// <Psi|M|Psi>_p
	auto Mcon = TreeFunctions::contraction(Psi, Mrep, tree); /// <Psi|M|Psi>_(p)
	cout << "Hrep:" << endl;
	Mrep.print();
	cout << "Hcon:\n";
	Mcon.print();
}

int main() {
	tree_examples();
	paper_1();
	return EXIT_SUCCESS;
}

