//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "TensorTree.h"

template <typename T>
TensorTree<T>::TensorTree(const TTBasis& basis) {
	Initialize(basis);
}

template <typename T>
void TensorTree<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node& node : basis) {
		attributes.emplace_back(Tensor<T>(node.TDim()));
	}
}

template <typename T>
void TensorTree<T>::Write(ostream& os) const {

	os.write("TTre", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const Tensor<T>& Phi : attributes) {
		os << Phi;
	}
	os << flush;
}

template <typename T>
void TensorTree<T>::Read(istream& is) {

	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TTre");
	assert(s_check == s_key);

	// Read number of nodes
	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes.clear();
	for (int i = 0; i < nnodes; i++) {
		Tensor<T> Phi(is);
		attributes.emplace_back(Phi);
	}
}

template <typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t) {
	t.Write(os);
	return os;
}

template <typename T>
istream& operator>>(istream& is, TensorTree<T>& t) {
	t.Read(is);
	return is;
}
