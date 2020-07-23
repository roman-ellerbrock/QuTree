//
// Created by Roman Ellerbrock on 2/10/20.
//

#ifndef TREEIO_H
#define TREEIO_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/SparseMatrixTree.h"

namespace TreeIO {
	void Status(size_t it, size_t max, size_t freq, size_t length);

	void StatusTime(size_t it, size_t max, size_t freq, size_t length,
		chrono::high_resolution_clock::time_point& t1, chrono::high_resolution_clock::time_point& t2,
		chrono::microseconds& time);

	void Output(const TensorTreecd& Psi, const Tree& tree, ostream& os = cout);

	template <typename T>
	void Occupancy(const TensorTree<T>& Psi, const Tree& tree, ostream& os = cout);

	template <typename T>
	void Leafs(const TensorTree<T>& Psi, const MatrixTree<T>& Rho, const Tree& tree, ostream& os = cout);

	template <typename T>
	Matrix<T> LeafDensity(const TensorTree<T>& Psi, const MatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree);

	template <typename T>
	Matrix<T> LeafDensity(const TensorTree<T>& Psi, const SparseMatrixTree<T>& Rho,
		const Leaf& leaf, const Tree& tree);

	template <class A>
	void print(const vector<A>& vec);
}

#endif //TREEIO_H
