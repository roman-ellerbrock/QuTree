#include "DenseOverlap.h"

template<typename T>
DenseOverlap<T>::DenseOverlap(istream& is) {
	Read(is);
}

template<typename T>
DenseOverlap<T>::DenseOverlap(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
DenseOverlap<T>::DenseOverlap(const TTBasis& basis) {
	Initialize(basis);
}

template<typename T>
DenseOverlap<T>::DenseOverlap(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {
	Initialize(basis);
	Calculate(Psi, Chi, basis);
}

template<typename T>
void DenseOverlap<T>::Initialize(const TTBasis& basis) {
	// Clear the overlaps for reinitialization
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		const size_t dim = tdim.getntensor();
		const auto k = (size_t) node.ChildIdx();
		FactorMatrix<T> mat(dim, k);
		attributes.push_back(mat);
	}
}

template<typename T>
FactorMatrix<T> DenseOverlap<T>::Calculate(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {

	for (const Node& node : basis) {
		CalculateLayer(Psi[node], Chi[node], node);
	}

	const Node& topnode = basis.TopNode();
	return this->operator[](topnode);
}

template<typename T>
void DenseOverlap<T>::CalculateLayer(const Tensor<T>& Phi,
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
Tensor<T> DenseOverlap<T>::TransformTensor(const Tensor<T>& Phi,
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
void DenseOverlap<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info(os);
		(*this).operator[](node).print(os);
	}
}

template<typename T>
void DenseOverlap<T>::print(ostream& os) const {
	for (const auto& x : *this) {
		x.print(os);
	}
}

template <typename T>
void DenseOverlap<T>::Write(ostream& os) const {
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
void DenseOverlap<T>::Read(istream& is) {
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
ostream& operator<<(ostream& os, const DenseOverlap<T>& S) {
	if (&os == &std::cout ) {
		S.print(os);
	} else{
		S.Write(os);
	}
	return os;
}

template <typename T>
istream& operator>>(istream& is, DenseOverlap<T>& S) {
	S.Read(is);
	return is;
}

typedef complex<double> cd;
template class DenseOverlap<cd>;
template ostream& operator<< <cd> (ostream&, const DenseOverlap<cd>& );
template istream& operator>> <cd> (istream&, DenseOverlap<cd>& );

typedef double d;
template class DenseOverlap<d>;
template ostream& operator<< <d> (ostream&, const DenseOverlap<d>& );
template istream& operator>> <d> (istream&, DenseOverlap<d>& );


