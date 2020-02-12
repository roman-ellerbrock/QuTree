//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "MatrixTree.h"

template<typename T>
MatrixTree<T>::MatrixTree(const TTBasis& basis) {
	Initialize(basis);
}

template <typename T>
void MatrixTree<T>::Initialize(const TTBasis& basis) {
	attributes_.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		attributes_.emplace_back(Matrix<T>(tdim.GetNumTensor()), tdim.GetNumTensor());
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

template<typename T>
void MatrixTree<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info();
		attributes_[node].print();
	}
}

template<typename T>
void MatrixTree<T>::print(ostream& os) const {
	for (const auto& mat : *this) {
		mat.print();
	}
}

