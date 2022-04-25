//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"
#include "Tensor/TensorLapack.h"

typedef complex<double> cd;

typedef double d;

template<typename T>
TensorTree<T>::TensorTree(const Tree& tree,
	function<Tensor<T>(const TensorShape&)> gen)
	:SubTreeAttribute<Tensor<T>>(tree) {

	for (const Node *node : this->nodes_) {
		(*this)[node] = gen(node->shape_);
	}

	for (const Edge* edge : this->edges_) {
		(*this)[edge] = (*this)[edge->from()];
	}

	normalize();
}

template<typename T>
Tensor<T> normalize(const Tensor<T>& phi, const Edge *edge, double eps) {
	size_t k = edge->outIdx();
	return ::normalize(phi, k, eps);
}

template Tensor<cd> normalize(const Tensor<cd>& phi, const Edge *edge, double eps);
template Tensor<d> normalize(const Tensor<d>& phi, const Edge *edge, double eps);

template<typename T>
void TensorTree<T>::normalize(double eps) {

	for (const Edge *edge: this->edges_) {
		Tensor<T>& Q = (*this)[edge];
		/// normalize Q and project on previous Q
		Tensor<T> preQ = Q;
		Q = ::normalize(Q, edge, eps);
		auto R = contraction(Q, preQ, edge->outIdx());

		/// transform next node
		Tensor<T>& B = (*this)[edge->to()];
		B = matrixTensor(R, B, edge->inIdx());

		/// transform next edges
		auto postedges = postEdges(edge);
		for (const Edge& postedge : postedges) {
			Tensor<T>& BQ = (*this)[postedge];
			BQ = matrixTensor(R, BQ, edge->inIdx());
		}
	}
}

template<typename T>
void TensorTree<T>::print() const {
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

template<typename T>
ostream& output(ostream& os, const TensorTree<T>& A) {
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

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& A);
