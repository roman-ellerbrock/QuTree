#include <gtest/gtest.h>
#include "Tensor/ManybodyTensor.h"
#include "Util/statistics.h"

TEST (ManybodyTensor, slice) {
    size_t n = 3;
//    Tensord A = randomd({n, n, n});
    Tensord A = aranged({n, n, n});
    ManybodyTensor M;
    TensorShape start({1, 1, 1});
    TensorShape end = A.shape_;
    auto S = slice(A, start, end);
    EXPECT_EQ(8, S.shape_.totalDimension());
    EXPECT_EQ(13, S[0]);
    EXPECT_EQ(14, S[1]);
    EXPECT_EQ(16, S[2]);
    EXPECT_EQ(17, S[3]);
    EXPECT_EQ(22, S[4]);
    EXPECT_EQ(23, S[5]);
    EXPECT_EQ(25, S[6]);
    EXPECT_EQ(26, S[7]);
}

TEST (ManybodyTensor, slice2) {
    size_t n = 3;
    Tensord A = aranged({n, n, n});
    ManybodyTensor M;
    TensorShape start({0, 0, 0});
    TensorShape end({3, 1, 1});
    auto S = slice(A, start, end);
    EXPECT_EQ(3, S.shape_.totalDimension());
    EXPECT_EQ(0, S[0]);
    EXPECT_EQ(1, S[1]);
    EXPECT_EQ(2, S[2]);
}

TEST (ManybodyTensor, slice3) {
    size_t n = 3;
    Tensord A = aranged({n, n, n});
    ManybodyTensor M;
    TensorShape start({0, 0, 0});
    TensorShape end({1, 3, 1});
    auto S = slice(A, start, end);
    EXPECT_EQ(3, S.shape_.totalDimension());
    EXPECT_EQ(0, S[0]);
    EXPECT_EQ(3, S[1]);
    EXPECT_EQ(6, S[2]);
}

TEST (ManybodyTensor, slice4) {
    size_t n = 3;
    Tensord A = aranged({n, n, n});
    ManybodyTensor M;
    TensorShape start({0, 0, 0});
    TensorShape end({1, 1, 3});
    auto S = slice(A, start, end);
    EXPECT_EQ(3, S.shape_.totalDimension());
    EXPECT_EQ(0, S[0]);
    EXPECT_EQ(9, S[1]);
    EXPECT_EQ(18, S[2]);
}

TEST (ManybodyTensor, recompose) {
    size_t n = 10;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);
    auto Ares = M.recompose();
    EXPECT_NEAR(0., residual(A, Ares), 1e-12);
}

TEST (ManybodyTensor, Core) {
    size_t n = 10;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);

    // safe reference data
    Tensord C = M.C_;
    auto U = M.U_;

    M.slicing({4, 4, 4});

    /// Check core tensor
    auto Cref = slice(C, {6, 6, 6}, {10, 10, 10});
    EXPECT_EQ(TensorShape({4, 4, 4}), M.C_.shape_);
    EXPECT_NEAR(0., residual(M.C_, Cref), 1e-12);

    /// Check trafos
    for (size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(TensorShape({10, 4}), M.U_[k].shape_);
        auto Uref = slice(U[k], {0, 6}, {10, 10});
        EXPECT_NEAR(0., residual(Uref, M.U_[k]), 1e-12);
    }
}

TEST (ManybodyTensor, FirstOrderSlices) {
    size_t n = 10;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);

    // safe reference data
    Tensord C = M.C_;
    auto U = M.U_;

    M.slicing({4, 4, 4});

    /// Check core tensor slices
    auto Cref0 = slice(C, {0, 6, 6}, {6, 10, 10});
    EXPECT_EQ(TensorShape({6, 4, 4}), M.C1_[0].shape_);
    EXPECT_NEAR(0., residual(M.C1_[0], Cref0), 1e-12);

    auto Cref1 = slice(C, {6, 0, 6}, {10, 6, 10});
    EXPECT_EQ(TensorShape({4, 6, 4}), M.C1_[1].shape_);
    EXPECT_NEAR(0., residual(M.C1_[1], Cref1), 1e-12);

    auto Cref2 = slice(C, {6, 6, 0}, {10, 10, 6});
    EXPECT_EQ(TensorShape({4, 4, 6}), M.C1_[2].shape_);
    EXPECT_NEAR(0., residual(M.C1_[2], Cref2), 1e-12);

    /// Check trafos
    for (size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(TensorShape({10, 6}), M.U1_[k].shape_);
        auto Uref = slice(U[k], {0, 0}, {10, 6});
        EXPECT_NEAR(0., residual(Uref, M.U1_[k]), 1e-12);
    }
}

/*TEST (ManybodyTensor, contraction0) {
    Tensord A = randomd({10, 10, 10});
    ManybodyTensor M;
    M.hosvd(A);
    M.slicing({5, 5, 5});
    A = M.recompose();
    Matrixd Sref = contraction(A, A, 0);
    Sref.print();
    Matrixd S({10, 10});
    contraction(S, M, M, 0);

    EXPECT_NEAR(0., residual(S, Sref), 1e-12);
}*/

TEST (ManybodyTensor, recompose_after_slicing_order0) {
    size_t n = 10;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);
    auto idx = indexMapping(0, A.shape_);
    for (size_t I = 0; I < A.shape_.totalDimension(); ++I) {
        indexMapping(idx, I, A.shape_);
        if (idx[0] < 6 || idx[1] < 6 || idx[2] < 6) {
            M.C_[I] = 0.;
        }
    }
    Tensord Aref = M.C_;
    for (size_t k = 0; k < M.U_.size(); ++k) {
        Aref = matrixTensor(M.U_[k], Aref, k);
    }
    M.slicing({4, 4, 4});
    Tensord Ares = M.C_;
    for (size_t k = 0; k < M.U_.size(); ++k) {
        Ares = matrixTensor(M.U_[k], Ares, k);
    }
    EXPECT_NEAR(0., residual(Aref, Ares), 1e-12);
}

TEST (ManybodyTensor, recompose_after_slicing_order1) {
    size_t n = 10;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);
    auto idx = indexMapping(0, A.shape_);
    for (size_t I = 0; I < A.shape_.totalDimension(); ++I) {
        indexMapping(idx, I, A.shape_);
        if (!(idx[0] < 6 && idx[1] >= 6 && idx[2] >= 6)) {
            M.C_[I] = 0.;
        }
    }

    Tensord Aref = M.C_;
    for (size_t k = 0; k < M.U_.size(); ++k) {
        Aref = matrixTensor(M.U_[k], Aref, k);
    }

    M.slicing({4, 4, 4});
    auto Ares = M.recompose();
    EXPECT_NEAR(0., residual(Aref, Ares), 1e-12);
}

TEST (ManybodyTensor, contraction) {
    size_t n = 10;
    size_t chi = 4;
    size_t delta = n - chi;
    Tensord A = randomd({n, n, n});

    ManybodyTensor M;
    M.hosvd(A);
    auto idx = indexMapping(0, A.shape_);
    for (size_t I = 0; I < A.shape_.totalDimension(); ++I) {
        indexMapping(idx, I, A.shape_);
        if (idx[0] < delta || idx[1] < delta || idx[2] < delta) {
            M.C_[I] = 0.;
        }
    }

    Tensord ref;
    {
        auto rho = contraction(M.C_, M.C_, 0);
        rho = M.U_[0] * rho * adjoint(M.U_[0]);
        auto tr = trace(rho);
        rho /= tr;
        auto rhox = heev(rho);
        ref = rhox.ev();
    }

    M.slicing({chi, chi, chi}, 0);
    Tensord res;
    {
        Matrixd rho({n, n});
        contraction(rho, M, M, 0);
        auto tr = trace(rho);
        rho /= tr;
        auto rhox = heev(rho);
        res = rhox.ev();
    }
    EXPECT_NEAR(0., residual(res, ref), 1e-12);
}

TEST (ManybodyTensor, contraction1) {
    size_t n = 10;
    size_t chi = 4;
    size_t delta = n - chi;
    Tensord A = randomd({n, n, n});
    ManybodyTensor M;
    M.hosvd(A);
    auto idx = indexMapping(0, A.shape_);
    for (size_t I = 0; I < A.shape_.totalDimension(); ++I) {
        indexMapping(idx, I, A.shape_);
        if (!(idx[0] < delta && idx[1] >= delta && idx[2] >= delta)) {
            M.C_[I] = 0.;
        }
    }

    Tensord ref;
    auto A2 = M.recompose();

    M.slicing({chi, chi, chi}, 1);

    {
        Tensord ref0 = contraction(A2, A2, 0);
        Tensord res0({n,n});
        contraction(res0, M, M, 0);
        EXPECT_NEAR(0., residual(res0, ref0), 1e-10);
    }

    {
        Tensord ref1 = contraction(A2, A2, 1);
        Tensord res1({n,n});
        contraction(res1, M, M, 1);
        EXPECT_NEAR(0., residual(res1, ref1), 1e-10);
    }
}

/*TEST (ManybodyTensor, spectrum) {
    size_t n = 400;
    for (size_t chi = 10; chi < 150; chi += 10) {
        size_t delta = n - chi;
        size_t n_samples = 3;
        vector<double> stat, stat0, stat1;
        size_t size, size0, size1;
        for (size_t i = 0; i < n_samples; ++i) {
            Tensord A = randomd({n, n, n});
            auto nrm = nrm2(A);
            A /= nrm;

            auto rho = contraction(A, A, 0);

            ManybodyTensor M0;
            M0.hosvd(A);
            M0.slicing({chi, chi, chi}, 0);

            Tensord rho0({n,n});
            contraction(rho0, M0, M0, 0);

            ManybodyTensor M1;
            M1.hosvd(A);
            M1.slicing({chi, chi, chi}, 1);

            Tensord rho1({n,n});
            contraction(rho1, M1, M1, 0);

            size = A.shape_.totalDimension();
            size0 = M0.size();
            size1 = M1.size();
            double avg = trace(rho);///(double) size;
            double avg0 = trace(rho0);///(double) M0.size();
            double avg1 = trace(rho1);///(double) M1.size();
            stat.push_back(avg);
            stat0.push_back(avg0);
            stat1.push_back(avg1);
        }

        auto eff = statistic_helper(stat);
        auto eff0 = statistic_helper(stat0);
        auto eff1 = statistic_helper(stat1);
        cout << n << "\t" << chi << "\t" << size0 << "\t" << eff0.first << "\t" << eff0.second;
        cout                     << "\t" << size1 << "\t" << eff1.first << "\t" << eff1.second << "\t" << eff.first << "\t" << eff.second << endl;
    }

    getchar();
}*/
