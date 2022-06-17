//
// Created by Roman Ellerbrock on 2/24/20.
//

#include <gtest/gtest.h>
#include "Util/GradientDescent_Implementation.h"
#include "Core/Vector.h"

class Interface {
public:

    Interface() {
        x_ = Vectord(2);
        x_(0) = 0.5;
        x_(1) = -0.2;
        alpha_ = 0.5;
        beta_ = 0.4;
    }

    ~Interface() = default;

    double evaluate(const Vectord& x) const {
        return -exp(-alpha_ * x(0) * x(0) - beta_ * x(1) * x(1)) + 1.;
    }

    Vectord& parameters() {
        return x_;
    }

    Vectord gradient() const {
        Vectord grad(2);
        grad(0) = 2 * alpha_ * x_(0) * evaluate(x_);
        grad(1) = 2 * beta_ * x_(1) * evaluate(x_);
        return grad;
    }

    void print() const {
        cout << "x = ";
        x_.print();
        cout << "f(x) = " << evaluate(x_) << endl;
    }

private:
    Vectord x_;
    double alpha_;
    double beta_;
};


TEST (GradientDescent, GradientDescent1) {
    Interface myInterface;
    gradientDescent<Interface, Vectord>(myInterface, 1., 20000);
    Vectord x = myInterface.parameters();
    ASSERT_NEAR(0., abs(x(0)), 0.05);
    ASSERT_NEAR(0., abs(x(1)), 0.05);
}
