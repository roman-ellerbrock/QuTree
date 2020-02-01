#include "FactorMatrixTree.h"

template<typename T>
FactorMatrixTree<T>::FactorMatrixTree(istream& is) {
	Read(is);
}

template<typename T>
FactorMatrixTree<T>::FactorMatrixTree(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
FactorMatrixTree<T>::FactorMatrixTree(const TTBasis& basis) {
	Initialize(basis);
}

template<typename T>
FactorMatrixTree<T>::FactorMatrixTree(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {
	Initialize(basis);
	Calculate(Psi, Chi, basis);
}

template<typename T>
void FactorMatrixTree<T>::Initialize(const TTBasis& basis) {
	// Clear the overlaps for reinitialization
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		const size_t dim = tdim.GetNumTensor();
		const auto k = (size_t) node.ChildIdx();
		FactorMatrix<T> mat(dim, k);
		attributes.push_back(mat);
	}
}

template<typename T>
FactorMatrix<T> FactorMatrixTree<T>::Calculate(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {

	for (const Node& node : basis) {
		CalculateLayer(Psi[node], Chi[node], node);
	}

	const Node& topnode = basis.TopNode();
	return this->operator[](topnode);
}

template<typename T>
void FactorMatrixTree<T>::CalculateLayer(const Tensor<T>& Phi,
	Tensor<T> AChi, const Node& node) {
	// Get references to the ACoefficients at each node
	if (!node.IsBottomlayer()) {
		for (int k = 0; k < node.nChildren(); k++) {
			// Get overlap-matrix from down_ under
			const Node& child = node.Down(k);
			const FactorMatrix<T>& spo = this->operator[](child);
			// Apply it to the right-hand side
			AChi = spo * AChi;
		}
	}
	int kchild = node.ChildIdx();
	Matrix<T> resultmat = Phi.DotProduct(AChi);
	this->operator[](node) = FactorMatrix<T>(resultmat, kchild);
}

template<typename T>
Tensor<T> FactorMatrixTree<T>::TransformTensor(const Tensor<T>& Phi,
	const Node& node) const {
	if (!node.IsBottomlayer()) {
		Tensor<T> hPhi(Phi);
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.Down(k);
			const FactorMatrix<T>& s = this->operator[](child);
			hPhi = s * hPhi;
		}
		return hPhi;
	} else {
		return Phi;
	}
}

/// I/O

template<typename T>
void FactorMatrixTree<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info(os);
		(*this).operator[](node).print(os);
	}
}

template<typename T>
void FactorMatrixTree<T>::print(ostream& os) const {
	for (const auto& x : *this) {
		x.print(os);
	}
}

template <typename T>
void FactorMatrixTree<T>::Write(ostream& os) const {
	os.write("TTDo", 4);

	// write number of nodes
	auto nnodes = (int32_t) attributes.size();
	os.write((char *) &nnodes, sizeof(int32_t));

	// Write Tensors
	for (const auto& Phi : attributes) {
		os << Phi;
	}
	os << flush;
}

template <typename T>
void FactorMatrixTree<T>::Read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TTDo");
	assert(s_key == s_check);

	int32_t nnodes;
	is.read((char *) &nnodes, sizeof(nnodes));

	// Read all Tensors
	attributes.clear();
	for (int i = 0; i < nnodes; i++) {
		FactorMatrix<T> M(is);
		attributes.emplace_back(M);
	}
}

template <typename T>
ostream& operator<<(ostream& os, const FactorMatrixTree<T>& S) {
	if (&os == &std::cout ) {
		S.print(os);
	} else{
		S.Write(os);
	}
	return os;
}

template <typename T>
istream& operator>>(istream& is, FactorMatrixTree<T>& S) {
	S.Read(is);
	return is;
}

typedef complex<double> cd;
template class FactorMatrixTree<cd>;
template ostream& operator<< <cd> (ostream&, const FactorMatrixTree<cd>& );
template istream& operator>> <cd> (istream&, FactorMatrixTree<cd>& );

typedef double d;
template class FactorMatrixTree<d>;
template ostream& operator<< <d> (ostream&, const FactorMatrixTree<d>& );
template istream& operator>> <d> (istream&, FactorMatrixTree<d>& );


