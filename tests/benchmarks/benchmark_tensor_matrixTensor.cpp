//
// Created by Roman Ellerbrock on 2/3/20.
//
#include <iomanip>
#include "benchmark_tensor_matrixTensor.h"
#include "benchmark_helper.h"

namespace benchmark {

	void matrixTensor_gemm(Tensorcd& C, const Matrixcd& h, const Tensorcd& B, size_t k) {
		if (k == 0) {

			size_t active = B.shape_[0];
			size_t activeC = C.shape_[0];
			size_t after = B.shape_.after(0);
			size_t m = activeC;
			size_t k = active; //activeB
			size_t n = after;

			blas::Layout layout = blas::Layout::ColMajor;
			blas::Op op_B = blas::Op::NoTrans;
			blas::Op op_h = blas::Op::Trans;
			blas::gemm(layout, op_h, op_B, m, n, k, 1., h.data(), m, B.data(), k,
				0., C.data(), m);
		} else if (k == B.shape_.lastIdx()) {

			size_t d = B.shape_.lastIdx();
			size_t m = C.shape_.before(d);
			size_t n = C.shape_[d];
			size_t k = B.shape_[d];

			blas::Layout layout = blas::Layout::ColMajor;
			blas::Op op_B = blas::Op::NoTrans;
			blas::Op op_h = blas::Op::Trans;
			blas::gemm(layout, op_B, op_h, 
				m, n, k, 
				1., 
				B.data(), m, 
				h.data(), n,
				0., 
				C.data(), m);
		} else {
			
			size_t before = B.shape_.before(k);
			size_t activeB = B.shape_[k];
			size_t activeC = C.shape_[k];
			size_t after = B.shape_.after(k);

			size_t preB = before * activeB;
			size_t preC = before * activeC;

			size_t m = before;
			size_t n = activeC;
			size_t k = activeB;

			blas::Layout layout = blas::Layout::ColMajor;
			blas::Op op_B = blas::Op::NoTrans;
			#pragma omp parallel for
			for (size_t aft = 0; aft < after; ++aft) {
				blas::gemm(layout, op_B, blas::Op::Trans, m, n, k, 
					1., 
					&B[aft * preB], m, 
					h.data(), n,
					0.,
					&C[aft * preC], m);
			}
		}
	}

	auto matrix_tensor_gemm(mt19937& gen, size_t dim, size_t order, size_t mode, 
	size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::high_resolution_clock> allocate_start, allocate_end;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Matrixcd S = randomcd({dim, dim});
		Tensorcd B = randomcd(tdim);
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			matrixTensor_gemm(B, S, A, mode);
			end = std::chrono::high_resolution_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto matrix_tensor(mt19937& gen, size_t dim, size_t order, size_t mode, 
	size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::high_resolution_clock> allocate_start, allocate_end;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A = randomcd(tdim);
		Matrixcd S = randomcd({dim, dim});
		Tensorcd B = randomcd(tdim);
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			matrixTensor(B, S, A, mode);
			end = std::chrono::high_resolution_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}


	void screen_dim_matrixTensor(mt19937& gen, ostream& os, size_t nsample, size_t mode, bool gemm_only, 
	size_t min_dim, size_t max_dim, size_t step_dim) {

		/// Screen dim of tensor
		size_t order = 3;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim, dim});
			gflops.push_back(calc_gflops<complex<double>>(shape, mode));
			os << dim << "\t" << gflops.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		cout << "# GB\n";
		vector<double> gB;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim, dim});
			gB.push_back(calc_mem_matrixTensor<complex<double>>(shape, mode));
			os << dim << "\t" << gB.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		if (!gemm_only) {
			// matrix Tensor
			size_t flopsidx = 0;
			os << "# matrix tensor product mode(" << mode << ") x" << nsample << "\n";
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		} else {
			// gemm only
			size_t flopsidx = 0;
			os << "# matrix tensor product mode(" << mode << ") - gemm only x" << nsample << "\n";
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = matrix_tensor_gemm(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9 << endl;
			}
		}

	}

	auto sample_peak_gemm(mt19937& gen, size_t dim, size_t nsample, ostream& os) {
		/// Initialize memory
		std::chrono::time_point<std::chrono::high_resolution_clock> allocate_start, allocate_end;
		Matrixcd A = randomcd({dim, dim});
		Matrixcd B = randomcd({dim, dim});
		Matrixcd C = randomcd({dim, dim});
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			gemm(C, A, B);
			end = std::chrono::high_resolution_clock::now();
			if (n == 0) { continue; }
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	void screen_order_matrixTensor(mt19937& gen, ostream& os, size_t nsample,
	 size_t which, bool gemm_only, 
     size_t mink, size_t maxk, size_t stepk) {
		/// Screen order of tensor
		size_t dim = 2;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t order = mink; order <= maxk; order += stepk) {
			size_t mode = order / 2 + 1;
			if (which == 0) {
				mode = 0;
			} else if (which == 2) {
				mode == order - 1;
			}
			TensorShape shape = make_TensorDim(order, dim);
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
			os << "# matrix tensor product x" << nsample << "\n";
			for (size_t order = mink; order <= maxk; order += stepk) {
				size_t mode = order / 2 + 1;
				if (which == 0) {
					mode = 0;
				} else if (which == 2) {
					mode == order - 1;
				}
				TensorShape shape = make_TensorDim(dim, order);
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
				os << gb[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		} else {
			size_t flopsidx = 0;
			os << "# matrix tensor product - gemm only x" << nsample << endl;
			for (size_t order = mink; order <= maxk; order += stepk) {
				size_t mode = order / 2 + 1;
				if (which == 0) {
					mode = 0;
				} else if (which == 2) {
					mode == order - 1;
				}
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = matrix_tensor_gemm(gen, dim, order, mode, nsample, cout);
				os << gb[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
			}
		}
	}

	void benchmark_peak_gemm(mt19937& gen, ostream& os, size_t nsample) {
		size_t min_dim = 7232;
		size_t max_dim = 12000;
		size_t step_dim = 256;
		size_t mode = 0;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim});
			gflops.push_back(calc_gflops<complex<double>>(shape, mode));
			os << dim << "\t" << gflops.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		cout << "# GB\n";
		vector<double> gB;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			TensorShape shape({dim, dim});
			gB.push_back(calc_mem_matrixTensor<complex<double>>(shape, mode));
			os << dim << "\t" << gB.back() << "\t" << 1e-9 * dim * dim * dim * dim * 6 << endl;
		}

		size_t flopsidx = 0;
		os << "# peak gemm x" << nsample << endl;
		for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
			os << std::setprecision(6);
			os << dim << "\t";
			auto stat = sample_peak_gemm(gen, dim, nsample, cout);
			os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-9 << "\t" << stat.second * 1e-9<< endl;
		}
	}
}

