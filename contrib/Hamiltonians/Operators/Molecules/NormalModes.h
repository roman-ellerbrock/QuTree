//
// Created by Roman Ellerbrock on 4/29/22.
//

#ifndef NORMALMODES_H
#define NORMALMODES_H
#include "TreeOperators/CoordinateTransformation.h"
#include "Core/Matrix.h"


class NormalModes : public CoordinateTransformation {
public:
	/**
	 * @param file_U Normal mode matrix (full NxN)
	 * @param file_x0 reference point int cartesian coordinates
	 * @param nmodes  number of active normal modes: U(NxN) -> U (N x nmodes)
	 */
	NormalModes(const string& file_U, const string& file_x0, size_t nmodes) {
		cout << "normal modes init.\n";
		ifstream is(file_U);
		size_t dim1 = 1;
		size_t dim2 = 1;
		is >> dim1 >> dim2;
		Matrixd U(dim1, dim2);
		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < dim2; ++j) {
				is >> U(j, i);
			}
		}
		U_ = Matrixd(dim1, nmodes);
		int shift = dim2 - nmodes;
		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < nmodes; ++j) {
				U_(i, j) = U(i, shift + j);
			}
		}

		/// read reference point
		ifstream is_x0(file_x0);
		x0_ = Vectord(dim1);
		for (size_t i = 0; i < dim1; ++i) {
			is_x0 >> x0_(i);
		}
	}

	Vectord transform(const Vectord& q) const override {
		Vectord x = U_ * q;
		x += x0_;
		return x;
	}

protected:
	Matrixd U_;
	Vectord x0_;
};


#endif //NORMALMODES_H
