//
// Created by Roman Ellerbrock on 6/30/22.
//

#include "PolymorphicMemory_CUDA.h"
#include <blas.hh>

namespace polymorphic {
//	blas::set_device( device );

	constexpr size_t device = 0;
	constexpr size_t batch_size = 1000;

	using d = double;
	using cd = complex<double>;

/*	template <> double *allocate<double, CUDA>(const size_t n) {
		return blas::device_malloc<double>(n);
	}

	template<>
	void copy<d, CUDA>(d *dest, const d *src, const size_t n) {
		blas::Queue queue( device, batch_size );
		blas::copy(n, src, 1, dest, 1, queue);
	}

	template<>
	void transfer<d, CPU, CUDA>(d *d_dest, const d* h_src, const size_t n, const size_t m) {
		size_t lda = n;
		blas::Queue queue( device, batch_size );
		blas::device_setmatrix(n, m, h_src, lda, d_dest, lda, queue);
	}

	template<>
	void transferFromDevice<d, CUDA, CPU>(d *h_dest, const d* d_src, const size_t n, const size_t m) {
		size_t lda = n;
		blas::Queue queue( device, batch_size );
		blas::device_getmatrix(n, m, d_src, lda, h_dest, lda, queue);
	}
*/
//	template class memory<complex<double>, CUDA>;
	template 
	class memory<double, CUDA>;
}
