//
// Created by Roman Ellerbrock on 4/4/22.
//

#include "ApplySCF.h"

typedef complex<double> cd;

typedef double d;

template<typename T>
double applyIteration(TensorTree<T>& HPsi, vector<TensorTree<T>>& mat,
	const TensorTree<T>& Psi, const SOP<T>& H) {

	double err = 0.;
	for (const Edge *edge : mat.front().edges_) {
		/// apply to edge
		const Node& node = edge->from();
		Tensor<T> HPsi_old = HPsi[node];
		HPsi[node] = Psi[node];

		apply(HPsi[node], mat, H, &node);

		/// check error
		auto x = contraction(HPsi_old, HPsi[node], edge->outIdx());
		double overlap = abs(nrm2(x));
		err += pow(abs(1. - overlap), 2);

		/// normalize
		HPsi[edge] = normalize(HPsi[node], edge, 1e-10);

		/// represent operator at node
		contraction(mat, HPsi[edge], Psi[edge], H, edge);
	}
	return err;
}

template double applyIteration(TensorTree<cd>& HPsi, vector<TensorTree<cd>>& Hmat,
	const TensorTree<cd>& Psi, const SOP<cd>& H);
template double applyIteration(TensorTree<d>& HPsi, vector<TensorTree<d>>& Hmat,
	const TensorTree<d>& Psi, const SOP<d>& H);

template<typename T>
double error(TensorTree<T> HPsi, vector<TensorTree<T>> mat,
	const TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree) {
	contraction(mat, HPsi, Psi, H);
	cout << "<Psi|CNot|Psi>=\n";
	for (const Node *node : HPsi.nodes_) {
		HPsi[node] = Psi[node];
		apply(HPsi[node], mat, H, node);
		auto x = contraction(HPsi[node], Psi[node]);
		x.print();
	}
	getchar();
	return 0.;
}

template double error(TensorTree<cd> HPsi, vector<TensorTree<cd>> mat,
	const TensorTree<cd>& Psi, const SOP<cd>& H, const Tree& tree);
template double error(TensorTree<d> HPsi, vector<TensorTree<d>> mat,
	const TensorTree<d>& Psi, const SOP<d>& H, const Tree& tree);

template<typename T>
void apply(TensorTree<T>& HPsi, vector<TensorTree<T>>& mat,
	const TensorTree<T>& Psi, const SOP<T>& H, size_t n_iter) {
	double eps = 1e-12;
	for (size_t i = 0; i < n_iter; ++i) {
		cout << "iteration = " << i + 1 << endl;
		double err = applyIteration(HPsi, mat, Psi, H);
		cout << "Avg. Error: " << err << endl;
		cout << endl;
		if (err < eps) { break; }
	}
}

template void apply(TensorTree<cd>& HPsi, vector<TensorTree<cd>>& mat,
	const TensorTree<cd>& Psi, const SOP<cd>& H, size_t n_iter);
template void apply(TensorTree<d>& HPsi, vector<TensorTree<d>>& mat,
	const TensorTree<d>& Psi, const SOP<d>& H, size_t n_iter);
