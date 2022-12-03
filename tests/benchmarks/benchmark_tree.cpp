//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_tree.h"
#include "Tree/Tree.h"
#include "Tree/TreeFactory.h"
#include "benchmark_helper.h"

namespace benchmark {

	auto sample_matrixtree(TensorTreecd& Rho, const TensorTreecd& Psi,
		const Tree& tree, size_t nsample) {

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			contraction(Rho, Psi, Psi);
//			Rho.Calculate(Psi, tree);
			end = std::chrono::system_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

	pair<double, double> benchmark_matrixtree(mt19937& gen,
	size_t dim, size_t nleaves, size_t nsample) {

		/// Initialize memory
		Tree tree = balancedTree(nleaves, dim, dim);
		TensorTreecd Psi(tree, randomcd);
		TensorTree Rho = matrixTreecd(tree);

		return sample_matrixtree(Rho, Psi, tree, nsample);
	}

	void screen_nleaves(mt19937& gen, ostream& os, size_t nsample,
	size_t min_order, size_t max_order, size_t stepk) {

		/// Screen dim of nleaves
		os << "# step_k : " << stepk << endl;
		os << "# Dense Tree Operation x" << nsample << endl;
		for (size_t l = min_order; l < max_order; l *= 2) {
			auto stat = benchmark_matrixtree(gen, 2, l, nsample);
			os << l << "\t\t" << stat.first * 1e-3 << "\t" 
				<< stat.second * 1e-3 << "\t" << 2 * l - 1 << endl;
		}
	}

		/*
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = benchmark::holematrixtree(gen,  dim, order, nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
		 */
/*
		os << "# Factor-matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::factormatrixtree(gen,  dim, order, nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
		*/
		/*
		os << "# Sparse factor matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::sparse_factormatrixtree(gen,  dim, order, 100 * nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
		*/

/*		os << "# Sparse hole matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::sparse_holematrixtree(gen, dim, order, 100 * nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}*/
//	}


/*	pair<double, double> factormatrixtree_sample(MatrixTreecd & fmat,
		const TensorTreecd& Psi, const Tree& tree, size_t nsample) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0;n<nsample;++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			dotProduct(fmat, Psi, Psi, tree);
//			fmat.Calculate(Psi, Psi, tree);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}

pair<double, double> factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree = TreeFactory::balancedTree(nleaves, dim, dim);
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
		TreeFunctions::represent(fmat, M, Psi, tree);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}
*/
/*
pair<double, double> sparse_factormatrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree = TreeFactory::balancedTree(nleaves, dim, dim);
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
		TreeFunctions::contraction(hole, Psi, fmat, tree);
		end = std::chrono::system_clock::now();
		duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
	}
	return statistic_helper(duration_vec);
}

pair<double, double> sparse_holematrixtree(mt19937& gen, size_t dim, size_t nleaves,
	size_t nsample, ostream& os) {
	/// Initialize memory
	Tree tree = TreeFactory::balancedTree(nleaves, dim, dim);
	TensorTreecd Psi(gen, tree);
	Matrixcd X(2, 2);
	X(0, 0) = 0.5;
	X(1, 1) = 0.5;
	LeafMatrixcd x(X);
	MLOcd M(x, 0);
	M.push_back(x, nleaves - 1);
	SparseMatrixTreecd fmat(M, tree);
	TreeFunctions::represent(fmat, M, Psi, tree);
	SparseMatrixTreecd hole(M, tree);
	return sparse_holematrixtree_sample(hole, fmat, M, Psi, tree, nsample);
}
*/
}
