#pragma once
#include "TreeOperators/Potential.h"

class NOClPotential
	: public Potential {
public:
	explicit NOClPotential(bool dissociation = true);
	~NOClPotential() override = default;
	double evaluate(const Vectord& Xv, size_t part) const override;

private:
	double evaluateGS(const Vectord& Xv, size_t part) const;
	double evaluateS1(const Vectord& Xv, size_t part) const;
	bool dissociation;
};

