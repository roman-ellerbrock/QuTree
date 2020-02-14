//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREE_H
#define MATRIXTREE_H
#include "NodeAttribute.h"
#include "TreeHandling/Tree.h"

template <typename T>
class MatrixTree: public NodeAttribute<Matrix<T>>{
	using NodeAttribute<Matrix<T>>::attributes_;
public:

	MatrixTree() = default;

	explicit MatrixTree(const Tree& tree);

	explicit MatrixTree(istream& is);

	explicit MatrixTree(const string& filename);

	~MatrixTree() = default;

	void Initialize(const Tree& tree);

	/// I/O
	/// Print human readable
	void print(const Tree& tree, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	/// Read in binary format
	void Read(istream& is);
	void Read(const string& filename);

	/// Write in binary format
	void Write(ostream& os) const;
	void Write(const string& filename) const;
};

template<typename T>
ostream& operator<<(ostream& os, const MatrixTree<T>& S);

template<typename T>
istream& operator>>(istream& is, MatrixTree<T>& S);

typedef MatrixTree<complex<double>> MatrixTreecd;

typedef MatrixTree<complex<double>> MatrixTreed;


#endif //MATRIXTREE_H
