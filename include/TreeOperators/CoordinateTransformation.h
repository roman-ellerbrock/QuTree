//
// Created by Roman Ellerbrock on 7/2/20.
//

#ifndef COORDINATETRANSFORMATION_H
#define COORDINATETRANSFORMATION_H
#include "Core/Vector.h"

class CoordinateTransformation {
public:
	CoordinateTransformation() = default;
	~CoordinateTransformation() = default;

	virtual Vectord Transform(const Vectord& x) const { return x; }
};


#endif //COORDINATETRANSFORMATION_H
