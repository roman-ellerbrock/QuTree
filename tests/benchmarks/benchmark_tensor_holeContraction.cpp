#include "benchmark_tensor_holeContraction.h"
#include "benchmark_helper.h"
#include <iomanip>

namespace benchmark {

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
			size_t A = ket.shape_.before(k);
			size_t BL = bra.shape_[k];
			size_t BR = ket.shape_[k];
			size_t C = ket.shape_.after(k);

			for (size_t c = 0; c < C; ++c) {
				blas::Layout layout = blas::Layout::ColMajor;
				blas::Op ct = blas::Op::ConjTrans;
				blas::Op notrans = blas::Op::NoTrans;
				blas::gemm(layout, ct, notrans, 
				BL, BR, A, 
				alpha, 
				&bra.data()[A * BL * c], A, 
				&ket.data()[A * BR * c], A, 
				1.,
				h.data(), BL);
			}
/*			size_t d = ket.shape_.lastIdx();
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
			h.data(), BL);*/
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
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			contraction_gemm(S, A, B, mode);
			end = std::chrono::high_resolution_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto hole_product(mt19937& gen, size_t dim, size_t order, size_t mode,
		size_t nsample, ostream& os) {

		/// Initialize memory
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Tensorcd B = randomcd(tdim);
		Matrixcd S({dim, dim});

		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			contraction(S, A, B, mode);
			end = std::chrono::high_resolution_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}
		return statistic_helper(duration_vec);
	}


	void screen_dim_contraction(mt19937& gen, ostream& os, 
	size_t nsample, size_t mode, bool gemm_only,
	size_t min_dim, size_t max_dim, size_t step_dim) {

		/// Screen dim of tensor
		size_t order = 3;

		vector<double> gflops;
		cout << "# GFLOPS\n";
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim, dim});
			gflops.push_back(calc_gflops<complex<double>>(shape, mode));
			os << dim << "\t" << gflops.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		cout << "# GB\n";
		vector<double> gB;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim, dim});
			gB.push_back(calc_mem<complex<double>>(shape, mode));
			os << dim << "\t" << gB.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}
		if (gemm_only) {
			size_t flopsidx = 0;
			os << "# tensor hole contraction - gemm only x" << nsample << endl;
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product_gemm(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		} else {
			size_t flopsidx = 0;
			os << "# tensor hole product x" << nsample << endl;
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		}
	}

	void screen_order_contraction(mt19937& gen, ostream& os, size_t nsample,
	 size_t which, bool gemm_only, size_t mink, size_t maxk, size_t stepk) {
		/// Screen order of tensor
		size_t dim = 2;
		size_t max_order = 29;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t order = mink; order <= maxk; order += stepk) {
			TensorShape shape = make_TensorDim(order, dim);
			size_t mode = order / 2 + 1;
			if (which == 0) {
				mode = 0;
			} else if (which == 2) {
				mode == order - 1;
			}
			gflops.push_back(calc_gflops<complex<double>>(shape, mode));
			os << order << "\t" << gflops.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		vector<double> gb;
		cout << "# Gb\n";
		for (size_t order = mink; order <= maxk; order += stepk) {
			size_t mode = order / 2 + 1;
			if (which == 0) {
				mode = 0;
			} else if (which == 2) {
				mode == order - 1;
			}
			TensorShape shape = make_TensorDim(order, dim);
			gb.push_back(calc_mem<complex<double>>(shape, mode));
			os << dim << "\t" << order << "\t" << gb.back() << endl;
		}

		if (!gemm_only) {
			size_t flopsidx = 0;
			os << "# contraction, which: " << which << " x" << nsample << endl;
			for (size_t order = mink; order <= maxk; order += stepk) {
				size_t mode = order / 2 + 1;
				if (which == 0) {
					mode = 0;
				} else if (which == 2) {
					mode == order - 1;
				}
				TensorShape shape = make_TensorDim(dim, order);
				os << std::setprecision(6);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product(gen, dim, order, mode, nsample, cout);
				os << gb[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		} else {
			size_t flopsidx = 0;
			os << "# Contraction - gemm only x" << nsample << ", which: " << which << endl;
			for (size_t order = mink; order <= maxk; order += stepk) {
				size_t mode = order / 2 + 1;
				if (which == 0) {
					mode = 0;
				} else if (which == 2) {
					mode == order - 1;
				}
				os << std::setprecision(6);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product_gemm(gen, dim, order, mode, nsample, cout);
				os << gb[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		}
	}
	
}