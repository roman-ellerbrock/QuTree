//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREE_H
#define MATRIXTREE_H
#include "NodeAttribute.h"
#include "TreeShape/Tree.h"

template <typename T>
class MatrixTree: public NodeAttribute<Matrix<T>>{
	using NodeAttribute<Matrix<T>>::attributes_;
public:

	MatrixTree() = default;

	explicit MatrixTree(const Tree& tree);

	MatrixTree(const Tree& tree1, const Tree& tree2);

	explicit MatrixTree(istream& is);

	explicit MatrixTree(const string& filename);

	~MatrixTree() = default;

	void Initialize(const Tree& tree);

	/// I/O
	/// Print human readable
	void print(const Tree& tree, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	/// read in binary format
	void read(istream& is);
	void read(const string& filename);

	/// Write in binary format
	void write(ostream& os) const;
	void write(const string& filename) const;
};

template<typename T>
ostream& operator<<(ostream& os, const MatrixTree<T>& S);

template<typename T>
istream& operator>>(istream& is, MatrixTree<T>& S);

typedef MatrixTree<complex<double>> MatrixTreecd;

typedef MatrixTree<double> MatrixTreed;


#endif //MATRIXTREE_H
