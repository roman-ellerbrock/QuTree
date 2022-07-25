#include "Tensor.h"
#include "TensorBLAS2.h"


template <typename T, typename U, template <typename> class Dev>
void cast(Tensor<T, Dev>& L, const Tensor<U, Dev>& R) {
    for (size_t I = 0; I < L.shape_.totalDimension(); ++I) {
        L(I) = (T) R(I);
    }
}

template <typename T, typename U, template <typename> class Dev>
void offDiagonal(Tensor<T, Dev>& off, const Tensor<U, Dev>& full) {
    off = full;
    for (size_t i = 0; i < off.shape_[0]; ++i) {
        off(i, i) = 0;
    }
}


template <typename T, typename U, template <typename> class Dev>
class mxpTensor {
    /**
     * @brief Mixed precision (mxp) Matrix for mxp gemm
     *
     * Note: 
     * - "lp" means low precision (e.g. float)
     * - "hp" means high precision (e.g. double)
     */
public:
    mxpTensor(const TensorShape& shape)
    : diag_({shape[0]}), off_(shape), offLP_(shape) {
    }

    mxpTensor(const Tensor<T, Dev>& A)
        : mxpTensor(A.shape_) {

        diag_ = diagonal(A);

        offDiagonal(off_, A);
    }

    void toLP() const {
        cast(offLP_, off_);
    }

    void addLP() {
        Tensor<T, Dev> highP(offLP_.shape_);
        cast(highP, offLP_);
        off_ += highP;
        offLP_.zero();
    }

    Tensor<T, Dev> convert() const {
        Tensor<T, Dev> full = off_;
        for (size_t i = 0; i < full.shape_[0]; ++i) {
            full(i, i) += diag_(i);
        }
        return full;
    }

    Tensor<T, Dev> diag_;
    Tensor<T, Dev> off_;
    mutable Tensor<U, Dev> offLP_;
};

using mxpTensord = mxpTensor<double, float, polymorphic::hostMemory>;

template <typename T, template <typename> class Dev>
void mdiagm(Tensor<T, Dev>& C, const Tensor<T, Dev>& B, const Tensor<T, Dev>& diag) {
    /**
     * @brief Multiply a dense matrix with a diagonal matrix
     * @param C output Matrix that is written on
     * @param diag diagonal Matrix
     * @param B dense input Matrix
     * 
     * C_ij = C_ij + B_ij * A_ii
     */
    if (B.shape_.order() != 2) {
        cerr << "Error: expected order of B to be two but received " << B.shape_.order() << "\n";
        exit(3);
    }
    if (diag.shape_.order() != 1) {
        cerr << "Error: expected order of diag to be one but received " << diag.shape_.order() << "\n";
        exit(3);
    }
    size_t n = diag.shape_[0];
    size_t inc_a = 1;
    size_t inc_b = 1;
    for (size_t i = 0; i < n; ++i) {
        T alpha = diag(i);
        const T* Bstart = B.data() + i * n;
        T* Cstart = C.data() + i * n;
	    blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b);
    }
}

template <typename T, template <typename> class Dev>
void diagmm(Tensor<T, Dev>& C, const Tensor<T, Dev>& diag, const Tensor<T, Dev>& B) {
    /**
     * @brief Multiply a dense matrix with a diagonal matrix
     * @param C output Matrix that is written on
     * @param B dense input Matrix
     * @param diag diagonal Matrix
     * 
     * C_ij = C_ij + A_ii * B_ij
     */
    if (B.shape_.order() != 2) {
        cerr << "Error: expected order of B to be two but received " << B.shape_.order() << "\n";
        exit(3);
    }
    if (diag.shape_.order() != 1) {
        cerr << "Error: expected order of diag to be one but received " << diag.shape_.order() << "\n";
        exit(3);
    }

    size_t n = diag.shape_[0];
    size_t inc_a = n;
    size_t inc_b = n;
    for (size_t i = 0; i < n; ++i) {
        T alpha = diag(i);
        const T* Bstart = B.data() + i;
        T* Cstart = C.data() + i;
	    blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b);
    }
}

template <typename T, typename U, template <typename> class Dev>
void gemm(mxpTensor<T, U, Dev>& C, const mxpTensor<T, U, Dev>& A,
const mxpTensor<T, U, Dev>& B) {
    /// n^1 diagonal * diagonal
    C.diag_ += productElementwise(A.diag_, B.diag_);

    /// n^2 diagonal * off & off & diagonal
    diagmm(C.off_, A.diag_, B.off_);
    mdiagm(C.off_, A.off_, B.diag_);

    /// perform n^3 in low precision
    A.toLP(); // off(hp) -> off(lp)
    B.toLP(); // off(hp) -> off(lp)
    gemm(C.offLP_, A.offLP_ , B.offLP_);
    C.addLP(); // off(hp) += off(lp)
}

template <typename T, typename U, template <typename> class Dev>
mxpTensor<T, U, Dev> gemm(const mxpTensor<T, U, Dev>& A, const mxpTensor<T, U, Dev>& B) {
    mxpTensor<T, U, Dev> C(A.off_.shape_);
    gemm(C, A, B);
    return C;
}
