//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "TensorTree.h"
#include "Core/Tensor_Extension.h"

template<typename T>
TensorTree<T>::TensorTree(const TTBasis& basis) {
	TensorTree::Initialize(basis);
}

template<typename T>
TensorTree<T>::TensorTree(istream& is) {
	Read(is);
}

template<typename T>
TensorTree<T>::TensorTree(const string& filename) {
	ifstream is(filename);
	Read(is);
}

/// Create tensor tree and occupy the coefficients
template<typename T>
TensorTree<T>::TensorTree(const TTBasis& basis,
	mt19937& gen, bool delta_lowest)
	: TensorTree(basis) {
	Generate(basis, gen, delta_lowest);
}

template<typename T>
void TensorTree<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node& node : basis) {
		attributes.emplace_back(Tensor<T>(node.TDim()));
	}
}

template<typename T>
void TensorTree<T>::Generate(const TTBasis& basis, mt19937& gen, bool delta_lowest) {
	assert(basis.nNodes() == attributes.size());
	for (const Node& node : basis) {
		if (node.IsBottomlayer()) {
			FillBottom(this->operator[](node), node);
		} else {
			FillUpper(this->operator[](node), gen, node, delta_lowest);
		}
	}
}

template<typename T>
void TensorTree<T>::FillUpper(Tensor<T>& Phi,
	mt19937& gen, const Node& node, bool delta_lowest) {

	assert(Phi.Dim().GetDimTot() > 0);
	Tensor_Extension::Generate(Phi, gen);
	// Set ground-state to "Hartree-Product" if flag is set
	if (delta_lowest) {
		for (size_t i = 0; i < Phi.Dim().GetDimPart(); ++i) {
			Phi(i, 0) = 0.;
		}
		Phi(0, 0) = 1.;
	}

	// orthonormalize
	GramSchmidt(Phi);
}

template<typename T>
void TensorTree<T>::FillBottom(Tensor<T>& Phi,
	const Node& node) {
	const Leaf& coord = node.PhysCoord();
	const PrimitiveBasis& grid = coord.PrimitiveGrid();
	grid.InitSPF(Phi);
}

/// (File) I/O
template<typename T>
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

template<typename T>
void TensorTree<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template<typename T>
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

template<typename T>
void TensorTree<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info(os);
		this->operator[](node).print(os);
	}
}

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t) {
	t.Write(os);
	return os;
}

template<typename T>
istream& operator>>(istream& is, TensorTree<T>& t) {
	t.Read(is);
	return is;
}

