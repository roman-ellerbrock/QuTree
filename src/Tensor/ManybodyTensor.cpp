#include "ManybodyTensor.h"
#include "TensorBLAS2.h"

void contraction(Matrixd& S, const ManybodyTensor& A, const ManybodyTensor& B, size_t k) {
    /// Core-Core   
    {
        Tensord C = B.C_;
        for (size_t l = 0; l < A.shape_.order(); ++l) {
            if (l == k) { continue; }
            Matrixd s = adjoint(A.U_[l]) * B.U_[l];
            C = matrixTensor(s, C, l);
        }
        auto rho = contraction(A.C_, C, k);

        S = B.U_[k] * rho * adjoint(A.U_[k]);
    }
 
   { /// <B^l|B^l>_(k) | k != l
        for (size_t l = 0; l < A.U1_.size(); ++l) {
            Tensord B1 = B.C1_[l];
/*            if (l != k) {
                Matrixd s = adjoint(A.U1_[l]) * B.U1_[l];
                B1 = matrixTensor(s, B1, k);
            }
            for (size_t m = 0; m < A.U_.size(); ++m) {
                if (m == l || m == k) { continue; }
                Matrixd s = adjoint(A.U_[l]) * B.U_[l];
                B1 = matrixTensor(s, B1, k);
            }*/
            auto rho = contraction(A.C1_[l], B1, k);
            if (l == k) {
                rho = B.U1_[k] * rho * adjoint(A.U1_[k]);
            } else {
                rho = B.U_[k] * rho * adjoint(A.U_[k]);
            }
            S += rho;
        }
    }
}