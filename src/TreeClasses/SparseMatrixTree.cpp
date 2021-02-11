//
// Created by Roman Ellerbrock on 2/12/20.
//

#include "TreeClasses/SparseMatrixTree.h"

template<typename T>
void SparseMatrixTree<T>::initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node *const node_ptr : sparseTree()) {
		const Node& node = *node_ptr;
		size_t dim = node.shape().lastDimension();
		attributes_.emplace_back(Matrix<T>(dim, dim));
	}
}

template<typename T>
void SparseMatrixTree<T>::print(ostream& os) const {
	for (const Node *node_ptr : sparseTree()) {
		const Node& node = *node_ptr;
		node.info(os);
		this->operator[](node).print();
	}
}

template<typename T>
void SparseMatrixTree<T>::write(ostream& os) const {
	os.write("SMTr", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes_.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& m : attributes_) {
		os << m;
	}
	os << flush;
}

template<typename T>
void SparseMatrixTree<T>::write(const string& filename) const {
	ofstream os(filename);
	write(os);
}

template<typename T>
void SparseMatrixTree<T>::read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("SMTr");
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

template<typename T>
void SparseMatrixTree<T>::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

template<typename T>
ostream& operator>>(ostream& os, const SparseMatrixTree<T>& hmat) {
	hmat.write(os);
}

template<typename T>
istream& operator<<(istream& is, SparseMatrixTree<T>& hmat) {
	hmat.read(is);
}

template class SparseMatrixTree<complex<double>>;

template class SparseMatrixTree<double>;

