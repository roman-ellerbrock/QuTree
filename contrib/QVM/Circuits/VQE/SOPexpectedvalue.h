//
// Created by Roman Ellerbrock on 8/20/21.
//

#ifndef SOPEXPECTEDVALUE_H
#define SOPEXPECTEDVALUE_H

Matrixcd expectation(const SOPcd& H,
	const TensorTreecd& Psi, const Tree& tree);

double expectationVal(const SOPcd& H,
	const TensorTreecd& Psi, const Tree& tree);

#endif //SOPEXPECTEDVALUE_H
