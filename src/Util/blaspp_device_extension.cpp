#include <cublas.h>
#include <limits>
#include "blaspp_device_extension.h"

namespace blas {
namespace device {
     void axpy(int64_t n,
          float alpha, 
          float const *x, int64_t incx, 
          float       *y, int64_t incy,
          blas::Queue& queue ) {
          cublasSaxpy( queue.handle(), n, &alpha, x, incx, y, incy );
     }

     void axpy(int64_t n,
          double alpha, 
          double const *x, int64_t incx, 
          double       *y, int64_t incy,
          blas::Queue& queue ) {
          cublasDaxpy( queue.handle(), n, &alpha, x, incx, y, incy );
     }

     void axpy(int64_t n,
          std::complex<float> alpha, 
          std::complex<float> const *x, int64_t incx, 
          std::complex<float>       *y, int64_t incy,
          blas::Queue& queue ) {
          cublasCaxpy( queue.handle(), n, 
          (const cuComplex*) &alpha, 
          (const cuComplex*) x, incx, 
          (      cuComplex*) y, incy );
     }

     void axpy(int64_t n, 
          std::complex<double> alpha, 
          std::complex<double> const *x, int64_t incx, 
          std::complex<double>       *y, int64_t incy,
          blas::Queue& queue ) {
          cublasZaxpy( queue.handle(), 
          n, 
          (const cuDoubleComplex*) &alpha, 
          (const cuDoubleComplex*) x, incx, 
          (      cuDoubleComplex*) y, incy );
     }
}
}

namespace blas {
/// @ingroup axpy
template <typename T>
void axpy(
    int64_t n,
    T alpha,
    T const *x, int64_t incx,
    T       *y, int64_t incy,
    Queue& queue )
{
    // check arguments
    blas_error_if( n < 0 );      // standard BLAS returns, doesn't fail
    blas_error_if( incx == 0 );  // standard BLAS doesn't detect inc[xy] == 0
    blas_error_if( incy == 0 );

    // check for overflow in native BLAS integer type, if smaller than int64_t
    if (sizeof(int64_t) > sizeof(device_blas_int)) {
        blas_error_if( n              > std::numeric_limits<device_blas_int>::max() );
        blas_error_if( std::abs(incx) > std::numeric_limits<device_blas_int>::max() );
        blas_error_if( std::abs(incy) > std::numeric_limits<device_blas_int>::max() );
    }

    device_blas_int n_    = (device_blas_int) n;
    device_blas_int incx_ = (device_blas_int) incx;
    device_blas_int incy_ = (device_blas_int) incy;

    blas::set_device( queue.device() );
    device::axpy( n_, alpha, x, incx_, y, incy_, queue );
}

template  void axpy(
    int64_t n,
    float alpha,
    float const *x, int64_t incx,
    float       *y, int64_t incy,
    Queue& queue );

template  void axpy(
    int64_t n,
    double alpha,
    double const *x, int64_t incx,
    double       *y, int64_t incy,
    Queue& queue );

template  void axpy(
    int64_t n,
    std::complex<float> alpha,
    std::complex<float> const *x, int64_t incx,
    std::complex<float>       *y, int64_t incy,
    Queue& queue );

template  void axpy(
    int64_t n,
    std::complex<double> alpha,
    std::complex<double> const *x, int64_t incx,
    std::complex<double>       *y, int64_t incy,
    Queue& queue );

}