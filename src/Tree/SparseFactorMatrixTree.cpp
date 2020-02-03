#include "SparseFactorMatrixTree.h"

vector<size_t> cast_to_vector_size_t(const vector<int>& a) {
	vector<size_t> b(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		b[i] = a[i];
	}
	return b;
}

template<typename T>
void SparseFactorMatrixTree<T>::Initialize(const TTBasis& basis) {
	attributes_.clear();
	for (const Node *const node_ptr : Active()) {
		const Node& node = *node_ptr;
		size_t dim = node.TDim().GetNumTensor();
		attributes_.emplace_back(FactorMatrix<T>(dim, node.ChildIdx()));
	}
}

template<typename T>
void SparseFactorMatrixTree<T>::Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
	const MPO<T>& M, const TTBasis& basis) {
	for (size_t n = 0; n < Active().size(); ++n) {
		const Node& node = Active().MCTDHNode(n);
		CalculateLayer(Bra[node], Ket[node], M, node);
	}
}

template<typename T>
FactorMatrix<T> SparseFactorMatrixTree<T>::CalculateUpper(const Tensor<T>& Bra, const Tensor<T>& Ket,
	const Node& node) {
	// @TODO: Optimize with switchbool trick
	// Swipe through children and apply active_ children's SPOs.
	Tensor<T> hKet(Ket);
	for (size_t l = 0; l < node.nChildren(); l++) {
		const Node& child = node.Down(l);
		if (!Active(child)) { continue; }
		const FactorMatrix<T>& h = operator[](child);
		hKet = h * hKet;
	}

	// calculate overlap
	size_t ChildIdx = node.ChildIdx();
	Matrix<T> resultmat = Bra.DotProduct(hKet);
	return FactorMatrix<T>(resultmat, ChildIdx);
}

template<typename T>
FactorMatrix<T> SparseFactorMatrixTree<T>::CalculateBottom(const Tensor<T>& Bra,
	const Tensor<T>& Ket, const MPO<T>& M, const Node& node, const Leaf& phys) {

	Tensor<T> MKet = M.ApplyBottomLayer(Ket, phys);
	int ChildIdx = node.ChildIdx();
	Matrix<T> resultmat = Bra.DotProduct(MKet);
	return FactorMatrix<T>(resultmat, ChildIdx);
}

template<typename T>
void SparseFactorMatrixTree<T>::CalculateLayer(const Tensor<T>& Bra,
	const Tensor<T>& Ket, const MPO<T>& M, const Node& node) {
	if (!Active(node)) { return; }

	if (node.IsBottomlayer()) {
		operator[](node) = CalculateBottom(Bra, Ket, M, node, node.PhysCoord());
	} else {
		operator[](node) = CalculateUpper(Bra, Ket, node);
	}
}

template<typename T>
Tensor<T> SparseFactorMatrixTree<T>::Apply(const Tensor<T>& Phi, const MPO<T>& M,
	const Node& node) const {
	if (!Active(node)) { return Phi; }
	if (node.IsBottomlayer()) {
		const Leaf& phys = node.PhysCoord();
		return M.ApplyBottomLayer(Phi, phys);
	} else {
		return ApplyUpper(Phi, node);
	}
}

template<typename T>
Tensor<T> SparseFactorMatrixTree<T>::ApplyUpper(Tensor<T> Phi, const Node& node) const {
	Tensor<T> hPhi(Phi.Dim());
	bool switchbool = true;
	for (size_t k = 0; k < node.nChildren(); ++k) {
		const Node& child = node.Down(k);
		if (!Active(child)) { continue; }
		const FactorMatrix<T>& h = operator[](child);
		if (switchbool) {
			multAB(hPhi, h, Phi, true);
		} else {
			multAB(Phi, h, hPhi, true);
		}
		switchbool = !switchbool;
	}
	if (switchbool) {
		return Phi;
	} else {
		return hPhi;
	}
}

template<typename T>
Tensor<T> SparseFactorMatrixTree<T>::ApplyHole(Tensor<T> Phi, const Node& hole_node) const {
	// @TODO: Optimize with switchbool trick
	assert(!hole_node.IsToplayer());
	const Node& parent = hole_node.Up();
	size_t drop = hole_node.ChildIdx();

	for (size_t k = 0; k < parent.nChildren(); ++k) {
		const Node& child = parent.Down(k);
		if ((child.ChildIdx() == drop) || (!Active(child))) { continue; }
		const FactorMatrix<T>& h = operator[](child);
		const TensorDim& tdim = Phi.Dim();
		Phi = h * Phi;
	}
	return Phi;
}

/// I/O functionality

template<typename T>
void SparseFactorMatrixTree<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node *const node_ptr : Active()) {
		const Node& node = *node_ptr;
		node.info();
		this->operator[](node).print();
	}
}

template<typename T>
void SparseFactorMatrixTree<T>::Write(ostream& os) const {
	os.write("FMTr", 4);

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
void SparseFactorMatrixTree<T>::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

template<typename T>
void SparseFactorMatrixTree<T>::Read(istream& is) {
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("FMTr");
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

template<typename T>
void SparseFactorMatrixTree<T>::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

template<typename T>
ostream& operator>>(ostream& os, const SparseFactorMatrixTree<T>& hmat) {
	hmat.Write(os);
}

template<typename T>
istream& operator<<(istream& is, SparseFactorMatrixTree<T>& hmat) {
	hmat.Read(is);
}

template
class SparseFactorMatrixTree<complex<double>>;

template
class SparseFactorMatrixTree<double>;
