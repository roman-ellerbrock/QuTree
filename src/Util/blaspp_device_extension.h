#pragma once
#include <blas.hh>

namespace blas {
template <typename T>
void axpy(
    int64_t n,
    T alpha,
    T const *x, int64_t incx,
    T       *y, int64_t incy,
    Queue& queue );
}