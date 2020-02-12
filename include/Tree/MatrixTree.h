//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREE_H
#define MATRIXTREE_H
#include "TreeStructuredObject.h"
#include "TensorTreeBasis/TensorTreeBasis.h"

template <typename T>
class MatrixTree: public TreeStructuredObject<Matrix<T>>{
	using TreeStructuredObject<Matrix<T>>::attributes_;

	MatrixTree() = default;

	MatrixTree(const TTBasis& basis);

	~MatrixTree() = default;

	void Initialize(const TTBasis& basis);

	/// I/O
	/// Print human readable
	void print(const TTBasis& basis, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	/// Read in binary format
	void Read(istream& is);

	/// Write in binary format
	void Write(ostream& os) const;
};

template<typename T>
ostream& operator<<(ostream& os, const MatrixTree<T>& S);

template<typename T>
istream& operator>>(istream& is, MatrixTree<T>& S);

typedef MatrixTree<complex<double>> MatrixTreecd;

typedef MatrixTree<complex<double>> MatrixTreed;


#endif //MATRIXTREE_H
