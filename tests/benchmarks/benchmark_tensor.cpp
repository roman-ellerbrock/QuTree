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

	template <typename T>
	double calc_gflops(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		return 1e-9 * (mul_ops * n + add_ops * n);
	}

	template <typename T>
	double calc_mem(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		/// read/write h, read A & B
		return 1e-9 * (2 * dim[k] * dim[k] + 2 * dim.totalDimension()) * sizeof(T);
	}

	template <typename T>
	double calc_mem_matrixTensor(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		/// B = h * A
		/// read/write B, read h & A
		return 1e-9 * (dim[k] * dim[k] + 3 * dim.totalDimension()) * sizeof(T);
	}

	auto hole_product_sample(Matrixcd& S, const Tensorcd& A,
		const Tensorcd& B, size_t nsample, size_t k) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample + 1; ++n) {
			std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			contraction(S, A, B, k);
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

	void screen_order_matrixTensor(mt19937& gen, ostream& os, size_t nsample,
	 size_t which, bool gemm_only) {
		/// Screen order of tensor
		size_t dim = 2;
		size_t max_order = 29;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t order = 3; order <= max_order; order += 2) {
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
		for (size_t order = 3; order <= max_order; order += 2) {
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
			for (size_t order = 3; order <= max_order; order += 2) {
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
			for (size_t order = 3; order <= max_order; order += 2) {
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

	void screen_order_contraction(mt19937& gen, ostream& os, size_t nsample,
	 size_t which, bool gemm_only) {
		/// Screen order of tensor
		size_t dim = 2;
		size_t max_order = 29;

		cout << "# GFLOPS\n";
		vector<double> gflops;
		for (size_t order = 3; order <= max_order; order += 2) {
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
		for (size_t order = 3; order <= max_order; order += 2) {
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
			for (size_t order = 3; order <= max_order; order += 2) {
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
				os << gb[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
			}
		} else {
			size_t flopsidx = 0;
			os << "# Contraction - gemm only x" << nsample << ", which: " << which << endl;
			for (size_t order = 3; order <= max_order; order += 2) {
				size_t mode = order / 2 + 1;
				if (which == 0) {
					mode = 0;
				} else if (which == 2) {
					mode == order - 1;
				}
				os << std::setprecision(6);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product_gemm(gen, dim, order, mode, nsample, cout);
				os << gb[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
			}
		}
	}
	
	void screen_dim_matrixTensor(mt19937& gen, ostream& os, size_t nsample, size_t mode, bool gemm_only) {

		/// Screen dim of tensor
		size_t order = 3;
		size_t min_dim = 10;
		size_t max_dim = 1000;
		size_t step_dim = 100;

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
			size_t flopsidx = 0;
			os << "# matrix tensor product mode(" << mode << ") x" << nsample << "\n";
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = matrix_tensor(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
			}
		} else {
			// add gemm part
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

	void screen_dim_contraction(mt19937& gen, ostream& os, size_t nsample, size_t mode, bool gemm_only) {

		/// Screen dim of tensor
		size_t order = 3;
		size_t min_dim = 64;
		size_t max_dim = 640;
		size_t step_dim = 64;

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
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
			}
		} else {
			size_t flopsidx = 0;
			os << "# tensor hole product x" << nsample << endl;
			for (size_t dim = min_dim; dim <= max_dim; dim += step_dim) {
				os << std::setprecision(12);
				os << dim << "\t" << order << "\t";
				auto stat = hole_product(gen, dim, order, mode, nsample, cout);
				os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
			}
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
			os << gflops[flopsidx] << "\t" << gB[flopsidx++] << "\t" << stat.first * 1e-6 << "\t" << stat.second * 1e-6<< endl;
		}
	}
}

