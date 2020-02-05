//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_tree.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "Tree/HoleMatrixTree.h"
#include "benchmark_helper.h"
#include "Tree/SparseHoleMatrixTree.h"


namespace benchmark {
	auto holematrixtree_sample(HoleMatrixTreecd& Rho, const TensorTreecd& Psi,
		const TTBasis& basis, size_t nsample) {

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			Rho.Calculate(Psi, basis);
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
		HoleMatrixTreecd Rho(basis);

		return holematrixtree_sample(Rho, Psi, basis, nsample);
	}

	pair<double, double> factormatrixtree_sample(FactorMatrixTreecd & fmat,
		const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0;n<nsample;++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			fmat.Calculate(Psi, Psi, basis);
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
	FactorMatrixTreecd fmat(basis);
	return factormatrixtree_sample(fmat, Psi, basis, nsample);
}

auto sparse_factormatrixtree_sample(SparseFactorMatrixTreecd& fmat, const MPOcd& M,
	const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		fmat.Calculate(Psi, M, basis);
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
	SPOMcd x(X);
	MPOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseFactorMatrixTreecd fmat(M, basis);
	return sparse_factormatrixtree_sample(fmat, M, Psi, basis, nsample);
}

auto sparse_holematrixtree_sample(SparseHoleMatrixTreecd& hole,
	const SparseFactorMatrixTreecd& fmat, const MPOcd& M,
	const TensorTreecd& Psi, const TTBasis& basis, size_t nsample) {

	vector<chrono::microseconds> duration_vec;
	for (size_t n = 0; n < nsample; ++n) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		hole.Calculate(Psi, fmat, basis);
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
	SPOMcd x(X);
	MPOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseFactorMatrixTreecd fmat(Psi, M, basis);
	SparseHoleMatrixTreecd hole(M, basis);
	return sparse_holematrixtree_sample(hole, fmat, M, Psi, basis, nsample);
}

}
