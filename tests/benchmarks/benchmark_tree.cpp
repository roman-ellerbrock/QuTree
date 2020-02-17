//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_tree.h"
#include "TreeShape/Tree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "benchmark_helper.h"
#include "SparseMatrixTreeFunctions.h"


namespace benchmark {
	using namespace MatrixTreeFunctions;
	auto holematrixtree_sample(MatrixTreecd& Rho, const TensorTreecd& Psi,
		const Tree& tree, size_t nsample) {

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			Contraction(Rho, Psi, tree, true);
//			Rho.Calculate(Psi, tree);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

	pair<double, double> holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os) {

		/// Initialize memory
		Tree tree(nleaves, dim, dim);
		TensorTreecd Psi(gen, tree);
		MatrixTreecd Rho(tree);

		return holematrixtree_sample(Rho, Psi, tree, nsample);
	}

	pair<double, double> factormatrixtree_sample(MatrixTreecd & fmat,
		const TensorTreecd& Psi, const Tree& tree, size_t nsample) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0;n<nsample;++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			DotProduct(fmat, Psi, Psi, tree);
//			fmat.Calculate(Psi, Psi, tree);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

pair<double, double> factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree(nleaves, dim, dim);
	TensorTreecd Psi(gen, tree);
	MatrixTreecd fmat(tree);
	return factormatrixtree_sample(fmat, Psi, tree, nsample);
}

auto sparse_factormatrixtree_sample(SparseMatrixTreecd& fmat, const MLOcd& M,
	const TensorTreecd& Psi, const Tree& tree, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		SparseMatrixTreeFunctions::Represent(fmat, M, Psi, tree);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}

pair<double, double> sparse_factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree(nleaves, dim, dim);
	TensorTreecd Psi(gen, tree);
	Matrixcd X(2, 2);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseMatrixTreecd fmat(M, tree);
	return sparse_factormatrixtree_sample(fmat, M, Psi, tree, nsample);
}

auto sparse_holematrixtree_sample(SparseMatrixTreecd& hole,
	const SparseMatrixTreecd& fmat, const MLOcd& M,
	const TensorTreecd& Psi, const Tree& tree, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		SparseMatrixTreeFunctions::Contraction(hole, Psi, fmat, tree);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}

pair<double, double> sparse_holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree(nleaves, dim, dim);
	TensorTreecd Psi(gen, tree);
	Matrixcd X(2, 2);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseMatrixTreecd fmat(M, tree);
	SparseMatrixTreeFunctions::Represent(fmat, M, Psi, tree);
	SparseMatrixTreecd hole(M, tree);
	return sparse_holematrixtree_sample(hole, fmat, M, Psi, tree, nsample);
}

}
