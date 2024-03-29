//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "TensorTree.h"
#include "Core/Tensor_Extension.h"

template<typename T>
TensorTree<T>::TensorTree(const Tree& tree) {
	TensorTree<T>::initialize(tree);
}

template<typename T>
TensorTree<T>::TensorTree(istream& is) {
	read(is);
}

template<typename T>
TensorTree<T>::TensorTree(const string& filename) {
	ifstream is(filename);
	read(is);
}

/// Create tensor tree and occupy the coefficients
template<typename T>
TensorTree<T>::TensorTree(std::mt19937& gen, const Tree& tree, bool delta_lowest)
	: TensorTree(tree) {
	fillRandom(gen, tree, delta_lowest);
}

template<typename T>
void TensorTree<T>::initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		attributes_.emplace_back(Tensor<T>(node.shape()));
	}
}

template<typename T>
void TensorTree<T>::fillRandom(std::mt19937& gen, const Tree& tree, bool delta_lowest) {
	assert(tree.nNodes() == attributes_.size());
	for (const Node& node : tree) {
		Tensor<T>& Phi = this->operator[](node);
		if (node.isBottomlayer()) {
			fillBottom(Phi, node);
		} else {
			fillUpper(Phi, gen, node, delta_lowest);
		}
	}
}

template<typename T>
void TensorTree<T>::fillUpper(Tensor<T>& Phi,
	mt19937& gen, const Node& node, bool delta_lowest) {

	assert(Phi.shape().totalDimension() > 0);
	Tensor_Extension::generate(Phi, gen);
	// Set ground-state to "Hartree-Product" if flag is set
	if (delta_lowest) {
		for (size_t i = 0; i < Phi.shape().lastBefore(); ++i) {
			Phi(i, 0) = 0.;
		}
		Phi(0, 0) = 1.;
	}

	// orthonormalize
	gramSchmidt(Phi);
}

template<typename T>
void TensorTree<T>::fillBottom(Tensor<T>& Phi,
	const Node& node) {
	const Leaf& coord = node.getLeaf();
	const LeafInterface& grid = coord.interface();
	grid.initSPF(Phi);
}

/// (File) I/O
template<typename T>
void TensorTree<T>::write(ostream& os) const {

	os.write("TTre", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes_.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const Tensor<T>& Phi : attributes_) {
		os << Phi;
	}
	os.flush();
}

template<typename T>
void TensorTree<T>::write(const string& filename, bool append) const {
	auto mod = ios::in;
	if (append) {
		mod = ios::in | ios::app;
	}
	ofstream os(filename, mod);
	write(os);
}

template<typename T>
void TensorTree<T>::read(istream& is) {

	if (is.bad()) {
		cerr << "Stream bad." << endl;
		exit(2);
	}
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TTre");
	if (s_check != s_key) {
		cerr << "Wrong keyword in wavefunction" << endl;
		cerr << "Expected: " << s_key << endl;
		cerr << "Actual: " << s_check << endl;
		exit(2);
	}

	// Read number of nodes
	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes_.clear();
	for (int i = 0; i < nnodes; i++) {
		Tensor<T> Phi(is);
		attributes_.emplace_back(Phi);
	}
}

template<typename T>
void TensorTree<T>::print(const Tree& tree, ostream& os) const {
	for (const Node& node : tree) {
		node.info(os);
		this->operator[](node).print(os);
	}
}

template <typename T>
void orthogonal(TensorTree<T>& Psi, const Tree& tree) {
	//Bottom-Up-Sweep
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			Tensor<T>& Phi = Psi[node];
			Matrix<T> S = Phi.dotProduct(Phi);
			auto spec = diagonalize(S);
			const auto& trafo = spec.first;
			const auto& eigenval = spec.second;

			for (size_t j = 0; j < S.dim1(); j++)
				assert(eigenval(j) >= -1e-12);

			Matrix<T> SW(S.dim1(), S.dim2());
			for (size_t j = 0; j < S.dim1(); j++) {
				for (size_t k = 0; k < S.dim1(); k++) {
					for (size_t l = 0; l < S.dim1(); l++) {
						SW(j, k) += trafo(k, l) * sqrt(max(eigenval(l), 0.))
							* conj(trafo(j, l));
					}
				}
			}

			const Node& parent = node.parent();
			Psi[parent] = matrixTensor(SW, Psi[parent], node.childIdx());
			gramSchmidt(Phi);
		}
	}
}

template <typename T>
void qrOrthogonal(TensorTree<T>& Psi, const Tree& tree) {

	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			Tensor<T>& Phi = Psi[node];
			Tensor<T> nPhi = qr(Phi);
			auto S = nPhi.dotProduct(Phi);
			const Node& parent = node.parent();
			Psi[parent] = matrixTensor(S, Psi[parent], node.childIdx());
			Psi[node] = nPhi;
		}
	}

}


template <typename T>
void orthonormal(TensorTree<T>& Psi, const Tree& tree) {
	for (const Node& node : tree) {
		gramSchmidt(Psi[node]);
	}
}

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t) {
	t.write(os);
	return os;
}

template<typename T>
istream& operator>>(istream& is, TensorTree<T>& t) {
	t.read(is);
	return is;
}

