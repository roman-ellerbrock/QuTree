#pragma once
#include "Core/stdafx.h"
#include <functional>
#include "Core/Tensor.h"

/**
 * \defgroup Util
 * \brief This group includes common utilites in QuTree.
 */

template<class Q, class T, typename U>
class BS_integrator
/**
 * \class BS_integrator
 * \ingroup Util
 * \brief This is a Bulirsch-Stoer Integrator.
 *
 * The design of this class is loosely related to the implementation
 * of a numerical recipes book.
 * */
{
public:
	BS_integrator() {}

	BS_integrator(int dim, T& initializer) {
		Initialize(dim, initializer);
	}

	~BS_integrator() {}

	void Initialize(int dim, T& initializer) {
		// set control parameters
		dim_ = dim;
		maxuse_ = 7;
		max_ = 11;
		shrink_ = 0.95;
		grow_ = 1.2;

		// calculate sequence_
		sequence_.resize(max_);
		sequence_[0] = 2;
		sequence_[1] = 4;
		sequence_[2] = 6;
		for (int i = 3; i < max_; i++)
			sequence_[i] = sequence_[i - 2] * 2;
//			sequence_[i] = 2 * (i + 1);

		// allocate work objects
		T a(initializer);
		for (int i = 0; i < 5; i++)
			yvec_.push_back(a);
		for (int i = 0; i < maxuse_; i++)
			ytab_.push_back(a);
		for (int i = 0; i < max_; i++)
			xtab_.push_back(0);
	}

	void Integrate(T& y, double& x, double xend, double& dx, double eps,
		function<void(Q&, double, T&, T&)> ddx, function<double(Q&, T&, T&)> err, Q& container) {
		while (xend - x > dx * 1e-6) {
			if (x + dx > xend) { dx = xend - x + 1E-10; }
//			cout << "\t" << "t_bs=" << x << " dt=" << dx << endl;
			y = bsstep(y, x, dx, eps, ddx, err, container);
		}
	}

	void clearmemory() {
		// yvec_
		for (int i = 0; i < yvec_.size(); i++) {
			T& ynow = yvec_[i];
			for (int j = 0; j < dim_; j++) {
				//				ynow(j) = 1E66;
				ynow(j) = 0;
			}
		}

		// ytab_
		for (int i = 0; i < ytab_.size(); i++) {
			T& ynow = ytab_[i];
			for (int j = 0; j < dim_; j++) {
				//				ynow(j) = 1E66;
				ynow(j) = 0;
			}
		}
		for (int i = 0; i < xtab_.size(); i++) {
			//			xtab_[i] = 1E66;
			xtab_[i] = 0;
		}
	}

	T bsstep(T y, double& x, double& dx, double eps,
		function<void(Q&, double, T&, T&)> ddx, function<double(Q&, T&, T&)> err, Q& container) {

		// Clear memory from last step (only for debug purpose)
		clearmemory();

		// save 
		T& y0 = yvec_[0];
		T& y3 = yvec_[3];
		y0 = y;
		y3 = y;

		// calculate derivative and save on [5]
		ddx(container, x, yvec_[4], y);

		// Initialize 
		double h = dx;
		bool end = false;

		while (!end) {
			int i = 0;
			// interpolation scheme
//			cout << "Seq: " ;
			while (!(end || (i >= max_))) {
//				if (i > 0) { cout << ", "; }
				// Midpoint integration
//				cout << sequence_[i] << endl;
				mmid(yvec_, x, h / (1.0 * sequence_[i]), sequence_[i], ddx, container);

				// Extrapolation
				double x_estimate = pow(h / (1. * sequence_[i]), 2);
				rzextr(i, x_estimate, yvec_[0], yvec_[1], yvec_[2]);

				yvec_[0] = yvec_[1];
				// y1=y1+y2
				T& y1 = yvec_[1];
				T& y2 = yvec_[2];

				for (int j = 0; j < dim_; j++)
					y1(j) += y2(j);

				// Calculate error and adjust step size
				double error = err(container, yvec_[0], yvec_[1]);

				if (error < eps) {
					x += h;
					if (i == maxuse_ - 1) {
						dx = h * shrink_;
					} else {
						if (i == maxuse_ - 2) {
							dx = h * grow_;
						} else {
							dx = (h * sequence_[maxuse_ - 2]) / (1. * sequence_[i]);
						}
					}
					end = true;
				}
//				cout << "i = " << i <<  ", seq[i] = " << sequence_[i] << endl;
				i++;
			}

			// Decrease step size, if it is too large
			if (!end) {
				h = h / 4.;
				for (int i = 0; i < (max_ - maxuse_) / 2.; i++) {
					h /= 2.;
				}
				if (dx + h == dx) {
					cout << "BS integrator: step size underflow" << endl;
					getchar();
					exit(1);
					assert(0);
				}
			}
		}
		y = y0;
		return y;
	}

	void mmid(vector<T>& y, double xstart, double h, int step, function<void(Q&, double, T&, T&)> ddx,
		Q& container) {
		T& y0 = y[0];
		T& y1 = y[1];
		T& y2 = y[2];
		T& y3 = y[3];
		T& y4 = y[4];

		y[0] = y[3];
		for (int i = 0; i < dim_; i++) {
			y1(i) = y3(i) + h * y4(i);
		}

		double x = xstart + h;

		ddx(container, x, y[2], y[1]);

		for (int n = 1; n < step; n++) {
			for (int i = 0; i < dim_; i++) {
				y2(i) = y0(i) + 2 * h * y2(i);
			}
			y[0] = y[1];
			y[1] = y[2];

			x += h;

			ddx(container, x, y2, y1);
		}

		// y= (y+ y1 + h*y2)*0.5
		for (int i = 0; i < dim_; i++) {
			y0(i) = (y0(i) + y1(i) + h * y2(i)) / 2.;
		}
	}

	void rzextr(int i_estimate, double x_estimate, T& yest,
		T& ysav, T& dy) {

		int mi = 0;
		vector<double> fx;
		fx.resize(maxuse_);

		xtab_[i_estimate] = x_estimate;

		if (i_estimate == 0) {
			ysav = yest;
			dy = yest;
			ytab_[0] = yest;
		} else {
			if (i_estimate < maxuse_) {
				mi = i_estimate + 1;
			} else {
				mi = maxuse_;
			}

			for (int k = 1; k < mi; k++) {
				fx[k] = xtab_[i_estimate - k] / x_estimate;
			}


			for (int j = 0; j < dim_; j++) {
				U ddy;
				U yy = yest(j);
				T& ynow = ytab_[0];
				U v = ynow(j);
				U c = yy;
				ynow(j) = yy;
				for (int k = 1; k < mi; k++) {
					U bi = fx[k] * v;
					U b = bi - c;
					if (abs(b) != 0) {
						b = (c - v) / b;
						ddy = c * b;
						c = bi * b;
					} else {
						ddy = v;
					}
					T& ynex = ytab_[k];
					v = ynex(j);
					ynex(j) = ddy;
					yy = yy + ddy;
				}
				dy(j) = ddy;
				ysav(j) = yy;
			}
		}
	}

protected:

	vector<int> sequence_;
	vector<T> yvec_;
	vector<double> xtab_;
	vector<T> ytab_;
	int dim_;
	int max_, maxuse_;
	double shrink_, grow_;
};

