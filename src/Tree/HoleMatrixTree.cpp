#include "HoleMatrixTree.h"

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const MPO<T>& M, const TTBasis& basis, const string& filename)
	: HoleMatrixTree(M, basis) {
	Read(filename);
}

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const TensorTree<T>& Psi,
	const FactorMatrixTree<T>& hmat, const MPO<T>& M, const TTBasis& basis)
	: HoleMatrixTree(M, basis) {
	Calculate(Psi, Psi, hmat, basis);
}

template<typename T>
void HoleMatrixTree<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node *const node_ptr : Active()) {
		const Node& node = *node_ptr;
		size_t dim = node.TDim().GetNumTensor();
		attributes.emplace_back(Matrix<T>(dim, dim));
	}
}

// Calculate Hole-Matrices
template<typename T>
void HoleMatrixTree<T>::Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
	const FactorMatrixTree<T>& hmat, const TTBasis& basis) {

	// Swipe top-down_ but exclude topnode
	int sub_topnode = Active().size() - 2;
	for (int n = sub_topnode; n >= 0; --n) {
		const Node& node = Active().MCTDHNode(n);
		assert(Active(node));

		const Node& parent = node.Up();
		Tensor<T> hKet = hmat.ApplyHole(Ket[parent], node);
		if (!parent.IsToplayer()) {
			hKet = multStateAB(this->operator[](parent), hKet);
		}
		operator[](node) = HoleProduct(Bra[parent], hKet, node.ChildIdx());
	}
}

template<typename T>
Tensor<T> HoleMatrixTree<T>::Apply(const Tensor<T>& Phi, const Node& node) const {
	if (node.IsToplayer()) {
		return Phi;
	} else {
		assert(Active(node));
		return multStateAB(this->operator[](node), Phi);
	}
}

/// I/O

template<typename T>
void HoleMatrixTree<T>::print(ostream& os) {
	for (const Node *node_ptr : Active()) {
		const Node& node = *node_ptr;
		node.info(os);
		this->operator[](node).print();
	}
}

template<typename T>
void HoleMatrixTree<T>::Write(ostream& os) const {
	os.write("HMTr", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& m : attributes) {
		os << m;
	}
	os << flush;
}

template<typename T>
void HoleMatrixTree<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template<typename T>
void HoleMatrixTree<T>::Read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("HMTr");
	assert(s_key == s_check);

	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes.clear();
	for (int i = 0; i < nnodes; i++) {
		Matrix<T> M(is);
		attributes.emplace_back(M);
	}
}

template<typename T>
void HoleMatrixTree<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
ostream& operator>>(ostream& os, const HoleMatrixTree<T>& hmat) {
	hmat.Write(os);
}

template<typename T>
istream& operator<<(istream& is, HoleMatrixTree<T>& hmat) {
	hmat.Read(is);
}

template
class HoleMatrixTree<complex<double>>;

template
class HoleMatrixTree<double>;

