//
// Created by Roman Ellerbrock on 2/3/20.
//
#include <iomanip>
#include "benchmark_tensor.h"
#include "benchmark_helper.h"
#include <blas/flops.hh>

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

/*	auto vector_hole_product_sample(Matrixcd& S, const vector<Tensorcd>& A,
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
				contraction(Ss[k], sAs[k], B[k], bef, act, act, aft);
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

*/
	auto hole_product_sample(Matrixcd& S, const Tensorcd& A,
		const Tensorcd& B, size_t nsample, size_t k) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			contraction(S, A, B, k);
			end = std::chrono::system_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto hole_product(mt19937& gen, size_t dim, size_t order, size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Tensorcd B = randomcd(tdim);
		Matrixcd S({dim, dim});

		/// Run hole-product
		return hole_product_sample(S, A, B, nsample, mode);
	}

	template<typename T>
	void contraction_gemm(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
		size_t k, T alpha = (T) 1., T beta = (T) 0.) {
		
		if (k == 0) {
			size_t BR = ket.shape_[0];
			size_t BL = bra.shape_[0];
			size_t C = ket.shape_.after(0);
			blas::Layout layout = blas::Layout::RowMajor;
			blas::Op ct = blas::Op::ConjTrans;
			blas::Op no = blas::Op::NoTrans;
			Tensor<T> h_add({BR, BL});
			blas::gemm(layout, ct, no, 
			BL, BR, C, 
			alpha, 
			bra.data(), BL,
			ket.data(), BR, 
			beta,
			h.data(), BR);
		} else if (k == bra.shape_.lastIdx()) {
			size_t d = ket.shape_.lastIdx();
			size_t BR = ket.shape_[d];
			size_t BL = bra.shape_[d];
			size_t A = ket.shape_.before(d); // C = 1
			blas::Layout layout = blas::Layout::ColMajor;
			blas::Op ct = blas::Op::ConjTrans;
			blas::Op notrans = blas::Op::NoTrans;
			blas::gemm(layout, ct, notrans, 
			BL, BR, A, 
			alpha, 
			bra.data(), A, 
			ket.data(), A, 
			beta,
			h.data(), BL);
		} else {

			size_t A = ket.shape_.before(k);
			size_t B = ket.shape_[k];
			size_t B2 = bra.shape_[k];
			size_t C = ket.shape_.after(k);

			size_t AC = A * C;
			blas::Layout layout = blas::Layout::ColMajor;
			blas::Op ct = blas::Op::ConjTrans;
			blas::Op notrans = blas::Op::NoTrans;
			blas::gemm(layout, ct, notrans, 
			B2, B, AC, 
			alpha, 
			bra.data(), AC, 
			ket.data(), AC, 
			beta,
			h.data(), B2);
		}
	}

	auto hole_product_gemm(mt19937& gen, size_t dim, size_t order,
	 size_t mode, size_t nsample, ostream& os) {
		/// Initialize memory
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Tensorcd B = randomcd(tdim);
		Matrixcd S({dim, dim});

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			contraction_gemm(S, A, B, mode);
			end = std::chrono::system_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto matrix_tensor_sample(Tensorcd& B, const Matrixcd& S, const Tensorcd& A,
		size_t nsample, size_t k) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			matrixTensor(B, S, A, k);
			end = std::chrono::system_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto matrix_tensor(mt19937& gen, size_t dim, size_t order, size_t mode, 
	size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Tensorcd S = randomcd({dim, dim});
		Tensorcd B(tdim);

		return matrix_tensor_sample(B, S, A, nsample, mode);
	}
/*
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
	*/

	template <typename T>
	double calc_gflops(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		return 1e-9 * (mul_ops * n + add_ops * n);
	}

	void screen_dim(mt19937& gen, ostream& os, size_t nsample) {

		/// Screen dim of tensor
		size_t order = 3;
		size_t min_dim = 100;
		size_t max_dim = 600;
		size_t step_dim = 100;
		size_t mode = 1;

		vector<double> gflops;
		cout << "# GFLOPS\n";
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim, dim});
			gflops.push_back(calc_gflops<complex<double>>(shape, mode));
			os << dim << "\t" << gflops.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}
		size_t flopsidx = 0;

/*		os << "# matrix tensor product\n";
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			os << std::setprecision(6);
			os << dim << "\t" << order << "\t";
			auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
			os << gflops[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
		}*/

/*		flopsidx = 0;
		os << "# gemm x" << nsample << endl;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			os << std::setprecision(6);
			os << dim << "\t" << order << "\t";
			auto stat = hole_product_gemm(gen, dim, order, mode, nsample, cout);
			os << gflops[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
		}*/

		flopsidx = 0;
		os << "# tensor hole product x" << nsample << endl;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			os << std::setprecision(6);
			os << dim << "\t" << order << "\t";
			auto stat = hole_product(gen, dim, order, mode, nsample, cout);
			os << gflops[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
		}

	}

/*
	void screen_nleaves(mt19937& gen, ostream& os, size_t nsample) {

		/// Screen dim of nleaves
		os << "# hole-matrix tree\n";
		size_t dim = 2;
		auto max_order = (size_t) pow(2, 20);
		size_t min_order = 4;
		*/
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
}

