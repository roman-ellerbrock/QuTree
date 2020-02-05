//
// Created by Roman Ellerbrock on 2/3/20.
//
#include <iomanip>
#include "../tests/benchmarks/benchmark_tensor.h"
#include "../tests/benchmarks/benchmark_helper.h"
#include "benchmark_tree.h"

namespace benchmark {
	TensorDim make_TensorDim(size_t order, size_t dim) {
		assert(order > 1);
		assert(dim > 0);
		vector<size_t> dims;
		for (size_t i = 0; i < order - 1; ++i) {
			dims.emplace_back(dim);
		}
		return TensorDim(dims, dim);
	}

	auto hole_product_sample(Matrixcd& S, const Tensorcd& A,
		const Tensorcd& B, size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			TensorHoleProduct(S, A, B, bef, act, act, aft);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto hole_product(mt19937& gen, size_t dim, size_t order, size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A(tdim, false);
		Tensorcd B(tdim, false);
		Tensor_Extension::Generate(A, gen);
		Tensor_Extension::Generate(B, gen);
		Matrixcd S(dim, dim);
		size_t aft = tdim.After(mode);
		size_t act = tdim.Active(mode);
		size_t bef = tdim.Before(mode);

		/// Run hole-product
		return hole_product_sample(S, A, B, nsample, bef, act, aft);
	}

	auto matrix_tensor_sample(Tensorcd& B, const Matrixcd& S, const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			mattensor(B, S, A, bef, act, act, aft, true);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto matrix_tensor(mt19937& gen, size_t dim, size_t order, size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
		allocate_start = std::chrono::system_clock::now();
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A(tdim, false);
		Matrixcd S(dim, dim);
		Tensor_Extension::Generate(A, gen);
		Tensor_Extension::Generate(S, gen);
		Tensorcd B(tdim, true);
		size_t aft = tdim.After(mode);
		size_t act = tdim.Active(mode);
		size_t bef = tdim.Before(mode);

		return matrix_tensor_sample(B, S, A, nsample, bef, act, aft);
	}

	void screen_order(mt19937& gen, ostream& os, size_t nsample) {
		/// Screen order of tensor
		size_t dim = 2;
		size_t max_order = 29;
		os << "# Matrix tensor\n";
		for (size_t order = 3; order <= max_order; order += 2) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}

		os << "# tensor hole product\n";
		for (size_t order = 3; order <= max_order; order += 2) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = hole_product(gen, dim, order, mode, nsample, cout);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screen_dim(mt19937& gen, ostream& os, size_t nsample) {

		/// Screen dim of tensor
		size_t order = 3;
		size_t min_dim = 50;
		size_t max_dim = 250;
		size_t step_dim = 20;
		os << "# matrix tensor product\n";
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}

		os << "# tensor hole product\n";
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = hole_product(gen, dim, order, mode, nsample, cout);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screen_nleaves(mt19937& gen, ostream& os, size_t nsample) {

		/// Screen dim of nleaves
		os << "# hole-matrix tree\n";
		size_t dim = 2;
		auto max_order = (size_t) pow(2, 20);
		size_t min_order = 4;
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
		os << "# Sparse factor matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::sparse_factormatrixtree(gen,  dim, order, 100 * nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
		*/

		os << "# Sparse hole matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::sparse_holematrixtree(gen,  dim, order, 100 * nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void run() {
		mt19937 gen(1989);
		size_t nsample = 20;
		ostream& os = cout;

//		screen_order(gen, os, nsample);
//		screen_dim(gen, os, nsample);
		screen_nleaves(gen, os, nsample);
	}
}

