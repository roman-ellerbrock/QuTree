#include <gtest/gtest.h>
#include "Operator.h"

using namespace qutree;

TEST(Operator, LeafOperator) {
    Tensor x = tensorlib::rand({2, 2}, tensorlib::real_options());
    ProductOperator P;
    P[5] = x;

    ASSERT_NEAR(0., (P[5] - x).norm().item<double>(), 1e-14);
}

TEST(Operator, SOP) {
    Tensor x = tensorlib::rand({2, 2}, tensorlib::real_options());
    ProductOperator P;
    P[5] = x;
    SumOfProducts S;
    S.push_back({1., P});
    auto Q = S[0];
    ASSERT_NEAR(0., (Q.second[5] - x).norm().item<double>(), 1e-14);
    ASSERT_NEAR(1., real(Q.first), 1e-14);
}