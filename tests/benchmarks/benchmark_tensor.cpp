//
// Created by Roman Ellerbrock on 2/3/20.
//
#include <iomanip>
#include "benchmark_tensor.h"
#include "benchmark_helper.h"
#include "benchmark_tree.h"

namespace benchmark {
	TensorShape make_TensorDim(size_t order, size_t dim) {
		assert(order > 1);
		assert(dim > 0);
		vector<size_t> dims;
		for (size_t i = 0; i < order; ++i) {
			dims.emplace_back(dim);
		}
		return TensorShape(dims);
	}

	auto vector_hole_product_sample(Matrixcd& S, const vector<Tensorcd>& A,
		const vector<Tensorcd>& B, const vector<Matrixcd>& Su, size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<Tensorcd> sAs;
		vector<Matrixcd> Ss;
		for (size_t k = 0; k < A.size(); ++k) {
			sAs.push_back(A[k]);
			Ss.push_back(S);
		}
			vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			for (size_t k = 0; k < A.size(); ++k) {
				multStateArTB(sAs[k], Su[k], A[k]);
				TensorContraction(Ss[k], sAs[k], B[k], bef, act, act, aft);
			}
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto vector_hole_product(mt19937& gen, size_t dim, size_t order, size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		auto tdim = make_TensorDim(order, dim);
		vector<Tensorcd> As;
		vector<Tensorcd> Bs;
		vector<Matrixcd> Sus;
		size_t chunk = pow(2, 21);
		for (size_t k = 0; k < chunk; ++k) {
			Tensorcd A(tdim, false);
			Tensorcd B(tdim, false);
			Matrixcd Su(dim, dim);
			Tensor_Extension::generate(A, gen);
			Tensor_Extension::generate(B, gen);
			Tensor_Extension::generate(Su, gen);
			As.emplace_back(A);
			Bs.emplace_back(B);
			Sus.emplace_back(Su);
		}
		Matrixcd S(dim, dim);
		size_t aft = tdim.after(mode);
		size_t act = tdim[mode];
		size_t bef = tdim.before(mode);

		/// Run hole-product
		return vector_hole_product_sample(S, As, Bs, Sus, nsample, bef, act, aft);
	}

	auto hole_product_sample(Matrixcd& S, const Tensorcd& A,
		const Tensorcd& B, size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		size_t chunk = pow(2, 20);
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			for (size_t k = 0; k < chunk; ++k) {
				TensorContraction(S, A, B, bef, act, act, aft);
			}
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
		Tensor_Extension::generate(A, gen);
		Tensor_Extension::generate(B, gen);
		Matrixcd S(dim, dim);
		size_t aft = tdim.after(mode);
		size_t act = tdim[mode];
		size_t bef = tdim.before(mode);

		/// Run hole-product
		return hole_product_sample(S, A, B, nsample, bef, act, aft);
	}

	auto matrix_tensor_sample(Tensorcd& B, const Matrixcd& S, const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			matrixTensor(B, S, A, bef, act, act, aft, true);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto matrix_tensor(mt19937& gen, size_t dim, size_t order, size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A(tdim, false);
		Matrixcd S(dim, dim);
		Tensor_Extension::generate(A, gen);
		Tensor_Extension::generate(S, gen);
		Tensorcd B(tdim, true);
		size_t aft = tdim.after(mode);
		size_t act = tdim[mode];
		size_t bef = tdim.before(mode);

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

		max_order = 3;
		os << "# tree-simulated tensor hole product\n";
		for (size_t order = 3; order <= max_order; order += 2) {
			size_t mode = order / 2 + 1;
			os << std::setprecision(6);
			os << dim << "\t" << order;
			auto stat = vector_hole_product(gen, dim, order, mode, nsample, cout);
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

		os << "# Factor-matrix tree\n";
		for (size_t order = min_order; order <= max_order; order *= 2) {
			size_t mode = order / 2 + 1;
			os << dim << "\t" << order;
			auto stat = benchmark::factormatrixtree(gen,  dim, order, nsample, os);
			os << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
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

