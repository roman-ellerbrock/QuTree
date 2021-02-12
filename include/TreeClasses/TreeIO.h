//
// Created by Roman Ellerbrock on 2/10/20.
//

#ifndef TREEIO_H
#define TREEIO_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/SparseMatrixTree.h"

namespace TreeIO {
	void status(size_t it, size_t max, size_t freq, size_t length);

	void statusTime(size_t it, size_t max, size_t freq, size_t length,
		chrono::high_resolution_clock::time_point& t1, chrono::high_resolution_clock::time_point& t2,
		chrono::microseconds& time);

	void output(const TensorTreecd& Psi, const Tree& tree, ostream& os = cout);

	template <typename T>
	void occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os = cout);

	template <typename T>
	void leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os = cout);

	template <typename T>
	Matrix<T> leafDensity(const TensorTree<T>& Psi, const MatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree);

	template <typename T>
	Matrix<T> leafDensity(const TensorTree<T>& Psi, const SparseMatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree);

	template <typename T>
	void entropyMap(const TensorTree<T>& Psi, const Tree& tree);

	template <class A>
	void print(const vector<A>& vec);
}

#endif //TREEIO_H
