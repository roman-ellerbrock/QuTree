//
// Created by Roman Ellerbrock on 2/12/20.
//

#include "TreeClasses/SparseMatrixTree.h"

template<typename T>
void SparseMatrixTree<T>::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node *const node_ptr : Active()) {
		const Node& node = *node_ptr;
		size_t dim = node.TDim().LastActive();
		attributes_.emplace_back(Matrix<T>(dim, dim));
	}
}

template<typename T>
void SparseMatrixTree<T>::print(ostream& os) {
	for (const Node *node_ptr : Active()) {
		const Node& node = *node_ptr;
		node.info(os);
		this->operator[](node).print();
	}
}

template<typename T>
void SparseMatrixTree<T>::Write(ostream& os) const {
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
void SparseMatrixTree<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template<typename T>
void SparseMatrixTree<T>::Read(istream& is) {
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
void SparseMatrixTree<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
ostream& operator>>(ostream& os, const SparseMatrixTree<T>& hmat) {
	hmat.Write(os);
}

template<typename T>
istream& operator<<(istream& is, SparseMatrixTree<T>& hmat) {
	hmat.Read(is);
}

template class SparseMatrixTree<complex<double>>;

template class SparseMatrixTree<double>;

