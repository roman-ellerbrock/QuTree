//
// Created by Roman on 3/29/2019.
//
#include "HoleOverlap.h"


template<typename T>
HoleOverlap<T>::HoleOverlap(istream& is) {
	Read(is);
}

template<typename T>
HoleOverlap<T>::HoleOverlap(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
HoleOverlap<T>::HoleOverlap(const TTBasis& basis){
	Initialize(basis);
}

template<typename T>
HoleOverlap<T>::HoleOverlap(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const DenseOverlap<T>& S,
		const TTBasis& basis){
	Initialize(basis);
	Calculate(Psi1, Psi2, S, basis);
}

template<typename T>
void HoleOverlap<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		size_t n = tdim.GetNumTensor();
		attributes.emplace_back(Matrix<T>(n, n));
	}
}

template<typename T>
void HoleOverlap<T>::CalculateLayer(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const DenseOverlap<T>& S, const Node& node) {

	const Node& parent = node.Up();
	const Tensor<T>& Bra = Psi[parent];
	Tensor<T> Ket = Chi[parent];

	// Transform sublying basis
	auto child_idx = (size_t) node.ChildIdx();
	for (size_t k = 0; k < parent.nChildren(); ++k) {
		if (k != child_idx) {
			const Node& child = parent.Down(k);
			const FactorMatrix<T>& h = S[child];
			Ket = h * Ket;
		}
	}

	// Transform with upper hole matrix
	if (!parent.IsToplayer()) {
		const Matrix<T>& holemat = this->operator[](parent);
		Ket = multStateAB(holemat, Ket);
	}

	// Calculate hole-product and save to attributes
	this->operator[](node) = HoleProduct(Bra, Ket, child_idx);
}

template<typename T>
void HoleOverlap<T>::Calculate(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const DenseOverlap<T>& S,
		const TTBasis& basis) {
	// Calculate the holematrices for two different wavefunctions.
	// @TODO: This routines offers many ways for improvement. Go ahead!

	// Top-down_ swipe. Note that the top-node is excluded from swipe!
	for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
		const Node& node = basis.GetNode(n);
		CalculateLayer(Psi1, Psi2, S, node);
	}
}

template <typename T>
void HoleOverlap<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info(os);
		(*this)[node].print(os);
	}
}

template <typename T>
void HoleOverlap<T>::print(ostream& os) const {
	for (const auto& x : *this) {
		x.print(os);
	}
}

template <typename T>
void HoleOverlap<T>::Write(ostream& os) const {
	os.write("TTHP", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& x : attributes) {
		os << x;
	}
	os << flush;
}

template <typename T>
void HoleOverlap<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template <typename T>
void HoleOverlap<T>::Read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TTHP");
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

template <typename T>
void HoleOverlap<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template <typename T>
ostream& operator<<(ostream& os, const HoleOverlap<T>& H) {
	if (&os == &cout) {
		H.print(os);
	} else {
		H.Write(os);
	}
	return os;
}

template <typename T>
istream& operator>>(istream& is, HoleOverlap<T>& H) {
	H.Read(is);
	return is;
}

typedef complex<double> cd;
template class HoleOverlap<complex<double>>;
template ostream& operator<< <cd> (ostream&, const HoleOverlap<cd>& );
template istream& operator>> <cd> (istream&, HoleOverlap<cd>& );

typedef double d;
template class HoleOverlap<d>;
template ostream& operator<< <d> (ostream&, const HoleOverlap<d>& );
template istream& operator>> <d> (istream&, HoleOverlap<d>& );
