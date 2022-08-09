#ifndef TENSORMXP_H
#define TENSORMXP_H
#include "Tensor.h"
#include "TensorBLAS2.h"


template <typename T, typename U, template <typename> class Dev, class ...Queue>
class mxpTensor;

template <typename T, typename U>
using Tensorm = mxpTensor<T, U, polymorphic::hostMemory>;

template <typename T, typename U>
void gemm(Tensorm<T, U>& C, const Tensorm<T, U>& A, const Tensorm<T, U>& B);

template <typename T, typename U>
Tensorm<T, U> gemm(const Tensorm<T, U>& A, const Tensorm<T, U>& B);

template <typename T, typename U, template <typename> class Dev, class ...Queue>
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
    : diag_({shape[0]}), off_(shape), offLP_(shape), shape_(shape) {
    }

    mxpTensor(const Tensor<T, Dev>& A, Queue& ...queue)
        : mxpTensor(A.shape_) {

        diag_ = diagonal(A, queue...);

        offDiagonal(off_, A, queue...);

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

    Tensor<T, Dev> convert(Queue& ...queue) const {
        Tensor<T, Dev> full = off_;
        addDiagonal(full, diag_, (T) 1., queue...);
        return full;
    }

/*    mxpTensor operator*(const mxpTensor& R) {
        return gemm(*this, R);
    }*/

    Tensor<T, Dev> diag_;
    Tensor<T, Dev> off_;
    Tensor<U, Dev> offLP_;
    TensorShape shape_;
};

template <typename T, typename U>
Tensorm<T, U> operator*(const Tensorm<T, U>& L, const Tensorm<T, U>& R);

using Tensordf = mxpTensor<double, float, polymorphic::hostMemory>;

#endif // TENSORMXP_H