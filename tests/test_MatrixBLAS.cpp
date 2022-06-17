//
// Created by Roman Ellerbrock on 6/4/21.
//
#include "Core/MatrixBLAS.h"
#include <gtest/gtest.h>

TEST(MatrixBLAS, QR) {
    Matrixcd A(5, 3);
    for (size_t i = 0; i < 15; ++i) {
        A[i] = i;
    }
    Matrixcd Q(A);
    qrBLAS(Q, A);
    complex<double> x = 0.;
    for (size_t i = 0; i < 5; ++i) {
        x += conj(Q(i, 1)) * Q(i, 1);
    }
    ASSERT_DOUBLE_EQ(x.real(),1.);
    ASSERT_DOUBLE_EQ(x.imag(),0.);
}
