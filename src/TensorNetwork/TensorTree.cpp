//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"
#include "Tensor/TensorLapack.h"

typedef complex<double> cd;

typedef double d;

template<typename T, template <typename> class Mem>
TensorTree<T,Mem>::TensorTree(const Tree& tree,
	function<Tensor<T,Mem>(const TensorShape&)> gen)
	:SubTreeAttribute<Tensor<T,Mem>>(tree) {

	for (const Node *node : this->nodes_) {
		(*this)[node] = gen(node->shape_);
	}

	for (const Edge* edge : this->edges_) {
		(*this)[edge] = (*this)[edge->from()];
	}

	normalize();
}

template<typename T, template <typename> class Mem = polymorphic::hostMemory>
Tensor<T,Mem> normalize(const Tensor<T,Mem>& phi, const Edge *edge, double eps) {
	size_t k = edge->outIdx();
	return ::normalize(phi, k, eps);
}

template Tensor<cd> normalize(const Tensor<cd>& phi, const Edge *edge, double eps);
template Tensor<d> normalize(const Tensor<d>& phi, const Edge *edge, double eps);

template<typename T, template <typename> class Mem>
void TensorTree<T,Mem>::normalize(double eps) {

	for (const Edge *edge: this->edges_) {
		Tensor<T,Mem>& Q = (*this)[edge];
		/// normalize Q and project on previous Q
		Tensor<T,Mem> preQ = Q;
		Q = ::normalize(Q, edge, eps);
		auto R = contraction(Q, preQ, edge->outIdx());

		/// transform next node
		Tensor<T,Mem>& B = (*this)[edge->to()];
		B = matrixTensor(R, B, edge->inIdx());

		/// transform next edges
		auto postedges = postEdges(edge);
		for (const Edge& postedge : postedges) {
			Tensor<T,Mem>& BQ = (*this)[postedge];
			BQ = matrixTensor(R, BQ, edge->inIdx());
		}
	}
}

template<typename T, template <typename> class Mem>
void TensorTree<T,Mem>::print() const {
	cout << "Nodes:\n";
	for (const Node *node : SubTree::nodes_) {
		node->info();
		(*this)[node].print();
	}
	cout << "Edges:\n";
	for (const Edge *edge : SubTree::edges_) {
		edge->info();
		(*this)[edge].print();
	}
}

template
class TensorTree<double>;

template
class TensorTree<complex<double>>;

template<typename T, template <typename> class Mem>
ostream& output(ostream& os, const TensorTree<T,Mem>& A) {
	os << "Natural occupancies:\n";
	for (const Edge *edge : A.edges_) {
//		if (!edge->isUpEdge()) { continue; }
		const Node& node = edge->from();
		auto rho = contraction(A[node], A[node], edge->outIdx());
		auto x = heev(rho);
		edge->info();
		rho.print();
		for (size_t i = 0; i < x.ev().shape_.totalDimension(); ++i) {
			os << x.ev()[i] << " ";
		}
		os << endl;
	}
	return os;
}

template ostream& output(ostream& os, const TensorTree<cd>& A);
template ostream& output(ostream& os, const TensorTree<d>& A);

template<typename T, template <typename> class Mem>
ostream& operator<<(ostream& os, const TensorTree<T,Mem>& A);
