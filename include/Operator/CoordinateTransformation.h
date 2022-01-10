//
// Created by Roman Ellerbrock on 7/2/20.
//

#ifndef COORDINATETRANSFORMATION_H
#define COORDINATETRANSFORMATION_H
#include "Tensor/Tensor"

class CoordinateTransformation {
public:
	CoordinateTransformation() = default;
	~CoordinateTransformation() = default;

	virtual Vectord transform(const Vectord& x) const { return x; }
};


#endif //COORDINATETRANSFORMATION_H
