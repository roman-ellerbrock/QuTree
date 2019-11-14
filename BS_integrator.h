#pragma once
#include "stdafx.h"
#include <functional>
#include "Tensor.h"


/**
 * \class BS_integrator
 * \ingroup QD-lib
 * \brief This is a Bulirsch-Stoer Integrator.
 *
 * The design of this class is losely related to the implementation
 * of a numerical recipes book.
 *
 * */

template<class Q, class T, typename U>
class BS_integrator
{
public:
	BS_integrator(){}
	~BS_integrator(){}

	void Initialize(int dim_, double eps_, T& initializer)
	{
		// set control parameters
		dim = dim_;
		eps = eps_;
		maxuse = 7;
		max = 11;
		shrink = 0.95;
		grow = 1.2;

		// calculate sequence
		sequence.resize(max);
		sequence[0] = 2;
		sequence[1] = 4;
		sequence[2] = 6;
		for (int i = 3; i < max; i++)
			sequence[i] = sequence[i - 2] * 2;
//			sequence[i] = 2 * (i + 1);

		// allocate work objects
		T a(initializer);
		for (int i = 0; i < 5; i++)
			yvec.push_back(a);
		for (int i = 0; i < maxuse; i++)
			ytab.push_back(a);
		for (int i = 0; i < max; i++)
			xtab.push_back(0);

	}

	void Integrate(T& y, double& x, double xend, double& dx,
		function<void(Q&, double, T&, T&)> ddx, function<double(Q&, T&, T&)> err, Q& container)
	{
		while (xend - x > dx*1e-6)
		{
			if (x + dx > xend) { dx = xend - x + 1E-10; }
//			cout << "\t" << "t=" << x << " dt=" << dx << endl;
			y=bsstep(y, x, dx, ddx, err, container);
		}
	}

	void clearmemory()
	{
		// yvec
		for (int i = 0; i < yvec.size(); i++)
		{
			T& ynow = yvec[i];
			for (int j = 0; j < dim; j++)
			{
				//				ynow(j) = 1E66;
				ynow(j) = 0;
			}

		}

		// ytab
		for (int i = 0; i < ytab.size(); i++)
		{
			T& ynow = ytab[i];
			for (int j = 0; j < dim; j++)
			{
				//				ynow(j) = 1E66;
				ynow(j) = 0;
			}
		}
		for (int i = 0; i < xtab.size(); i++)
		{
			//			xtab[i] = 1E66;
			xtab[i] = 0;
		}
	}

	T bsstep(T y, double& x, double& dx,
		function<void(Q&, double, T&, T&)> ddx, function<double(Q&, T&, T&)> err, Q& container)
	{

		// Clear memory from last step (only for debug purpose)
		clearmemory();

		// save 
		T& y0=yvec[0];
		T& y3=yvec[3];
		y0 = y;
		y3 = y;

		// calculate derivative and save on [5]
		ddx(container, x, yvec[4], y);

		// Initialize 
		double h = dx;
		bool end = false;

		while (!end)
		{
			int i = 0;
			// interpolation scheme
			while (!(end || (i >= max)))
			{
				// Midpoint integration
				mmid(yvec, x, h / (1.0*sequence[i]), sequence[i], ddx, container);

				// Extrapolation
				double x_estimate = pow(h / (1.*sequence[i]), 2);
				rzextr(i, x_estimate, yvec[0], yvec[1], yvec[2]);

				yvec[0] = yvec[1];
				// y1=y1+y2
				T& y1 = yvec[1];
				T& y2 = yvec[2];
				#pragma omp for
				for (int j = 0; j < dim; j++)
					y1(j) += y2(j);

				// Calculate error and adjust step size
				double error = err(container, yvec[0], yvec[1]);

				if (error < eps)
				{
					x += h;
					if (i == maxuse - 1)
					{
						dx = h*shrink;
					}
					else
					{
						if (i == maxuse - 2)
						{
							dx = h*grow;
						}
						else
						{
							dx = (h*sequence[maxuse - 2]) / (1.*sequence[i]);
						}
					}
					end = true;
				}
				i++;
			}

			// Decrease step size, if it is too large
			if (!end)
			{
				h = h / 4.;
				for (int i = 0; i < (max - maxuse) / 2.; i++)
				{
					h /= 2.;
				}
				if (dx + h == dx)
				{
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
		Q& container)
	{
		T& y0 = y[0];
		T& y1 = y[1];
		T& y2 = y[2];
		T& y3 = y[3];
		T& y4 = y[4];

		y[0] = y[3];
		for (int i = 0; i < dim; i++)
		{
			y1(i) = y3(i) + h*y4(i);
		}

		double x = xstart + h;

		ddx(container, x, y[2], y[1]);

		for (int n = 1; n < step; n++)
		{
			for (int i = 0; i < dim; i++)
			{
				y2(i) = y0(i) + 2 * h*y2(i);
			}
			y[0] = y[1];
			y[1] = y[2];

			x += h;

			ddx(container, x, y2, y1);
		}

		// y= (y+ y1 + h*y2)*0.5
		for (int i = 0; i < dim; i++)
		{
			y0(i) = (y0(i) + y1(i) + h*y2(i)) / 2.;
		}

	}

	void rzextr(int i_estimate, double x_estimate, T& yest,
		T& ysav, T& dy)
	{

		int mi = 0;
		vector<double> fx;
		fx.resize(maxuse);

		xtab[i_estimate] = x_estimate;

		if (i_estimate == 0)
		{
			ysav = yest;
			dy = yest;
			ytab[0] = yest;
		}
		else
		{
			if (i_estimate < maxuse)
			{
				mi = i_estimate + 1;
			}
			else {
				mi = maxuse;
			}
			#pragma omp for
			for (int k = 1; k < mi; k++)
			{
				fx[k] = xtab[i_estimate - k] / x_estimate;
			}

			#pragma omp for
			for (int j = 0; j < dim; j++)
			{
				U ddy;
				U yy = yest(j);
				T& ynow = ytab[0];
				U v = ynow(j);
				U c = yy;
				ynow(j) = yy;
				for (int k = 1; k < mi; k++)
				{
					U bi = fx[k] * v;
					U b = bi - c;
					if (abs(b) != 0)
					{
						b = (c - v) / b;
						ddy = c*b;
						c = bi*b;
					}
					else
					{
						ddy = v;
					}
					T& ynex = ytab[k];
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

	vector<int> sequence;
	vector<T> yvec;
	vector<double> xtab;
	vector<T> ytab;
	int dim;
	int max, maxuse;
	double eps;
	double shrink, grow;
};

