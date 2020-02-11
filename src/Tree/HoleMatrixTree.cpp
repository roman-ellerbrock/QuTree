//
// Created by Roman on 3/29/2019.
//
#include "HoleMatrixTree.h"


template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(istream& is) {
	Read(is);
}

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const TTBasis& basis){
	Initialize(basis);
}

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const FactorMatrixTree<T>& S,
		const TTBasis& basis){
	Initialize(basis);
	Calculate(Psi1, Psi2, S, basis);
}

template<typename T>
HoleMatrixTree<T>::HoleMatrixTree(const TensorTree<T>& Psi, const TTBasis& basis) {
	Initialize(basis);
	Calculate(Psi, basis);
}

template<typename T>
void HoleMatrixTree<T>::Initialize(const TTBasis& basis) {
	attributes_.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		size_t n = tdim.GetNumTensor();
		size_t order = tdim.GetOrder();
		attributes_.emplace_back(FactorMatrix<T>(n, n, order));
	}
}

////////////////////////////////////////////////////////////////////////
/// Calculate for orthogonal tensor trees
////////////////////////////////////////////////////////////////////////

template<typename T>
void HoleMatrixTree<T>::CalculateLayer(const TensorTree<T>& Psi,
	const Node& node) {

	const Node& parent = node.Up();
	Tensor<T> Ket(parent.TDim());
	const Tensor<T>& Bra = Psi[parent];

	// Transform with upper hole matrix
	if (!parent.IsToplayer()) {
		const FactorMatrix<T>& holemat = this->operator[](parent);
		multStateAB(Ket, holemat, Bra);
	} else {
		Ket = Bra;
	}

	// Calculate hole-product and save to attributes_
	auto child_idx = (size_t) node.ChildIdx();
	FactorMatrix<T>& s = this->operator[](node);
	HoleProduct(s, Bra, Ket, child_idx);
	s.Mode() = node.TDim().GetOrder();
}

template<typename T>
void HoleMatrixTree<T>::Calculate(const TensorTree<T>& Psi1,
	const TTBasis& basis) {
	// Calculate the holematrices for two different wavefunctions.

	// Top-down_ swipe. Note that the top-node is excluded from swipe!
	for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
		const Node& node = basis.GetNode(n);
		CalculateLayer(Psi1, node);
	}
}

////////////////////////////////////////////////////////////////////////
/// Calculate for non-othorgonal tensor trees
////////////////////////////////////////////////////////////////////////

template<typename T>
void HoleMatrixTree<T>::CalculateLayer(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const FactorMatrixTree<T>& S, const Node& node) {

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
		const FactorMatrix<T>& holemat = this->operator[](parent);
		Ket = multStateAB(holemat, Ket);
//		Ket = holemat * Ket;
	}

	// Calculate hole-product and save to attributes_
	FactorMatrix<T>& s = this->operator[](node);
	HoleProduct(s, Bra, Ket, child_idx);
	s.Mode() = node.TDim().GetOrder();
}

template<typename T>
void HoleMatrixTree<T>::Calculate(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const FactorMatrixTree<T>& S,
		const TTBasis& basis) {
	// Calculate the holematrices for two different wavefunctions.

	// Top-down_ swipe. Note that the top-node is excluded from swipe!
	for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
		const Node& node = basis.GetNode(n);
		CalculateLayer(Psi1, Psi2, S, node);
	}
}

////////////////////////////////////////////////////////////////////////
/// I/O functions
////////////////////////////////////////////////////////////////////////

template <typename T>
void HoleMatrixTree<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info(os);
		(*this)[node].print(os);
	}
}

template <typename T>
void HoleMatrixTree<T>::print(ostream& os) const {
	for (const auto& x : *this) {
		x.print(os);
	}
}

template <typename T>
void HoleMatrixTree<T>::Write(ostream& os) const {
	os.write("TTHP", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes_.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& x : attributes_) {
		os << x;
	}
	os << flush;
}

template <typename T>
void HoleMatrixTree<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template <typename T>
void HoleMatrixTree<T>::Read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TTHP");
	assert(s_key == s_check);

	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes_.clear();
	for (int i = 0; i < nnodes; i++) {
		FactorMatrix<T> M(is);
		attributes_.emplace_back(M);
	}
}

template <typename T>
void HoleMatrixTree<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template <typename T>
ostream& operator<<(ostream& os, const HoleMatrixTree<T>& H) {
	if (&os == &cout) {
		H.print(os);
	} else {
		H.Write(os);
	}
	return os;
}

template <typename T>
istream& operator>>(istream& is, HoleMatrixTree<T>& H) {
	H.Read(is);
	return is;
}

typedef complex<double> cd;
template class HoleMatrixTree<complex<double>>;
template ostream& operator<< <cd> (ostream&, const HoleMatrixTree<cd>& );
template istream& operator>> <cd> (istream&, HoleMatrixTree<cd>& );

typedef double d;
template class HoleMatrixTree<d>;
template ostream& operator<< <d> (ostream&, const HoleMatrixTree<d>& );
template istream& operator>> <d> (istream&, HoleMatrixTree<d>& );
