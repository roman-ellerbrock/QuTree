#include <gtest/gtest.h>
#include "Util/blaspp_device_extension.h"
#include "Tensor/cuTensor.h"
#include "Tensor/TensorBLAS1.h"

using namespace polymorphic;

template <typename T>
class BLASpp : public ::testing::Test {
    public:
    BLASpp() : n_(100), hx_(arange<T, hostMemory>({100})), hy_(arange<T, hostMemory>({100})) {
        x_ = transfer<T, cuMemory, hostMemory>(hx_);
        y_ = transfer<T, cuMemory, hostMemory>(hy_);
    }

    Tensor<T> hx_;
    Tensor<T> hy_;

    cuTensor<T> x_;
    cuTensor<T> y_;
    size_t n_;
};

using MyTypes = ::testing::Types<float, double, complex<double>>;
TYPED_TEST_SUITE(BLASpp, MyTypes);

TYPED_TEST(BLASpp, cuDaxpy) {
    TypeParam alpha = 1.;
    int incx = 1;
    int incy = 1;
    blas::Queue queue(0, 1000);

    blas::axpy(this->n_, 
    alpha,
    this->x_.data(), incx,
    this->y_.data(), incy,
    queue );

    blas::axpy(this->n_, 
    alpha,
    this->hx_.data(), incx,
    this->hy_.data(), incy );

    EXPECT_EQ(100, this->hy_.shape_.front());
    EXPECT_EQ(100, this->hx_.shape_.front());
    EXPECT_EQ(100, this->x_.shape_.front());
    EXPECT_EQ(100, this->y_.shape_.front());
    Tensor<TypeParam> hy2 = transfer<TypeParam, hostMemory, cuMemory>(this->y_);
    EXPECT_NEAR(0., residual(this->hy_, hy2), 1e-12);
}
