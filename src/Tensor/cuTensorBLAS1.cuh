#ifndef CUTENSORBLAS1_CUH
#define CUTENSORBLAS1_CUH
//#include <cublas.h>
#include <stdafx.h>

template <typename TL, typename TR>
__global__ void cudaCast(TL* dest, const TR* src, size_t n);

#endif // CUTENSORBLAS1_CUH