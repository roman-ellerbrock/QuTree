//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "TreeClasses/MatrixTree.h"

template<typename T>
MatrixTree<T>::MatrixTree(istream& is) {
	Read(is);
}

template<typename T>
MatrixTree<T>::MatrixTree(const string& filename) {
	Read(filename);
}

template<typename T>
MatrixTree<T>::MatrixTree(const Tree& tree) {
	Initialize(tree);
}

template <typename T>
void MatrixTree<T>::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorDim& tdim = node.TDim();
		size_t dim = tdim.LastActive();
		attributes_.emplace_back(Matrix<T>(dim, dim));
	}
}

template <typename T>
void MatrixTree<T>::Write(ostream& os) const {
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
void MatrixTree<T>::Write(const string& filename) const  {
	ofstream os(filename);
	Write(os);
}

template <typename T>
void MatrixTree<T>::Read(istream& is) {
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
void MatrixTree<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
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
	S.Write(os);
}

template<typename T>
istream& operator>>(istream& is, MatrixTree<T>& S) {
	S.Read(is);
}

template class MatrixTree<complex<double>>;
template class MatrixTree<double>;
