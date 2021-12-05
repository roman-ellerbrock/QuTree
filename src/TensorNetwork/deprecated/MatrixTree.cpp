//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "TreeClasses/MatrixTree.h"

template<typename T>
MatrixTree<T>::MatrixTree(istream& is) {
	read(is);
}

template<typename T>
MatrixTree<T>::MatrixTree(const string& filename) {
	read(filename);
}

template<typename T>
MatrixTree<T>::MatrixTree(const Tree& tree) {
	Initialize(tree);
}

template<typename T>
MatrixTree<T>::MatrixTree(const Tree& tree1, const Tree& tree2) {
	attributes_.clear();

	if (tree1.nNodes() != tree2.nNodes()) {
		cerr << "Error: in MatrixTree: trees do not share the same topology!\n";
		exit(1);
	}
	for (const Node& node1 : tree1) {
		const Node& node2 = tree2.getNode(node1.address());
		const TensorShape& tdim1 = node1.shape();
		const TensorShape& tdim2 = node2.shape();
		size_t dim1 = tdim1.lastDimension();
		size_t dim2 = tdim2.lastDimension();
		attributes_.emplace_back(Matrix<T>(dim1, dim2));
	}
}


template <typename T>
void MatrixTree<T>::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorShape& tdim = node.shape();
		size_t dim = tdim.lastDimension();
		attributes_.emplace_back(Matrix<T>(dim, dim));
	}
}

template <typename T>
void MatrixTree<T>::write(ostream& os) const {
	os.write("MTre", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes_.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& Phi : attributes_) {
		os << Phi;
	}
	os << flush;
}

template <typename T>
void MatrixTree<T>::write(const string& filename) const  {
	ofstream os(filename);
	write(os);
}

template <typename T>
void MatrixTree<T>::read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("MTre");
	assert(s_key == s_check);

	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes_.clear();
	for (int i = 0; i < nnodes; i++) {
		Matrix<T> M(is);
		attributes_.emplace_back(M);
	}
}

template <typename T>
void MatrixTree<T>::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

template<typename T>
void MatrixTree<T>::print(const Tree& tree, ostream& os) const {
	for (const Node& node : tree) {
		node.info();
		this->operator[](node).print();
	}
}

template<typename T>
void MatrixTree<T>::print(ostream& os) const {
	for (const auto& mat : *this) {
		mat.print();
	}
}

template<typename T>
ostream& operator<<(ostream& os, const MatrixTree<T>& S) {
	S.write(os);
}

template<typename T>
istream& operator>>(istream& is, MatrixTree<T>& S) {
	S.read(is);
}

template class MatrixTree<complex<double>>;
template class MatrixTree<double>;
