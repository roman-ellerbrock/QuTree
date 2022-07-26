#include "Tensor.h"
#include "TensorBLAS2.h"


template <typename T, typename U, template <typename> class Dev>
class mxpTensor;

template <typename T, typename U, template <typename> class Dev>
void gemm(mxpTensor<T, U, Dev>& C, const mxpTensor<T, U, Dev>& A,
const mxpTensor<T, U, Dev>& B);

template <typename T, typename U, template <typename> class Dev>
mxpTensor<T, U, Dev> gemm(const mxpTensor<T, U, Dev>& A, const mxpTensor<T, U, Dev>& B);

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

        toLP();
    }

    void toLP() {
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

    mxpTensor operator*(const mxpTensor& R) {
        return gemm(*this, R);
    }

    Tensor<T, Dev> diag_;
    Tensor<T, Dev> off_;
    Tensor<U, Dev> offLP_;
};

using mxpTensord = mxpTensor<double, float, polymorphic::hostMemory>;

/**
 * @bTODO: for cuTensordf
 * - cast f - d
 * - productElementwise
 */