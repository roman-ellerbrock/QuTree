//
// Created by Roman Ellerbrock on 4/15/20.
//

#ifndef QVMPARAMETERS_H
#define QVMPARAMETERS_H


struct QVMParameters {

	QVMParameters() : eps(1e-6), max_spf(50), plus_spf(4) {}
	double eps;
	size_t max_spf;
	size_t plus_spf;

};


#endif //QVMPARAMETERS_H
