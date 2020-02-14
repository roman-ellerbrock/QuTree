//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_tree.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "Tree/MatrixTreeFunctions.h"
#include "benchmark_helper.h"
#include "SparseMatrixTreeFunctions.h"


namespace benchmark {
	using namespace MatrixTreeFunctions;
	auto holematrixtree_sample(MatrixTreecd& Rho, const TensorTreecd& Psi,
		const TTBasis& basis, size_t nsample) {

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			Contraction(Rho, Psi, basis, true);
//			Rho.Calculate(Psi, basis);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

	pair<double, double> holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
		size_t nsample, ostream& os) {

		/// Initialize memory
		TTBasis basis(nleaves, dim, dim);
		TensorTreecd Psi(basis, gen);
		MatrixTreecd Rho(basis);

		return holematrixtree_sample(Rho, Psi, basis, nsample);
	}

	pair<double, double> factormatrixtree_sample(MatrixTreecd & fmat,
		const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0;n<nsample;++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			DotProduct(fmat, Psi, Psi, basis);
//			fmat.Calculate(Psi, Psi, basis);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

pair<double, double> factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	TTBasis basis(nleaves, dim, dim);
	TensorTreecd Psi(basis, gen);
	MatrixTreecd fmat(basis);
	return factormatrixtree_sample(fmat, Psi, basis, nsample);
}

auto sparse_factormatrixtree_sample(SparseMatrixTreecd& fmat, const MLOcd& M,
	const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		SparseMatrixTreeFunctions::Represent(fmat, M, Psi, basis);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}

pair<double, double> sparse_factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	TTBasis basis(nleaves, dim, dim);
	TensorTreecd Psi(basis, gen);
	FactorMatrixcd X(2, 1);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseMatrixTreecd fmat(M, basis);
	return sparse_factormatrixtree_sample(fmat, M, Psi, basis, nsample);
}

auto sparse_holematrixtree_sample(SparseMatrixTreecd& hole,
	const SparseMatrixTreecd& fmat, const MLOcd& M,
	const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		SparseMatrixTreeFunctions::Contraction(hole, Psi, fmat, basis);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}

pair<double, double> sparse_holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	TTBasis basis(nleaves, dim, dim);
	TensorTreecd Psi(basis, gen);
	FactorMatrixcd X(2, 1);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseMatrixTreecd fmat(M, basis);
	SparseMatrixTreeFunctions::Represent(fmat, M, Psi, basis);
	SparseMatrixTreecd hole(M, basis);
	return sparse_holematrixtree_sample(hole, fmat, M, Psi, basis, nsample);
}

}
