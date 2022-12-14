//
// Created by Roman Ellerbrock on 5/3/22.
//

#ifndef PORTFOLIOOPTIMIZATION_H
#define PORTFOLIOOPTIMIZATION_H
#include "TreeOperators/SumOfProductsOperator.h"

SOPcd meanVarianceAnalysis(const Tensord& mu, const Tensord& cov,
	size_t Na, size_t Nt, size_t Nq,
	double alpha, double gamma, double rho, double K);

	SOPcd meanVarianceAnalysis(string tickers, size_t Na, size_t Nt,
	size_t NaTot, size_t NtTot, size_t Nq,
	double alpha, double gamma, double rho, double K);

#endif //PORTFOLIOOPTIMIZATION_H

